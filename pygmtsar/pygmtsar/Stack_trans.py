# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_align import Stack_align
from .tqdm_dask import tqdm_dask

class Stack_trans(Stack_align):
    
    def define_trans_grid(self, coarsen):
        import numpy as np
        # select radar coordinates extent
        rng_max, yvalid, num_patch = self.PRM().get('num_rng_bins', 'num_valid_az', 'num_patches')
        azi_max = yvalid * num_patch
        #print ('azi_max', azi_max, 'rng_max', rng_max)
        # this grid covers the full interferogram area
        # common single pixel resolution
        #azis = np.arange(0, azi_max+coarsen[0], coarsen[0], dtype=np.int32)
        #rngs = np.arange(0, rng_max+coarsen[1], coarsen[1], dtype=np.int32)
        # for subpixel resolution
        #azis = np.arange(0, azi_max+coarsen[0], coarsen[0], dtype=np.float64)
        #rngs = np.arange(0, rng_max+coarsen[1], coarsen[1], dtype=np.float64)
        # this grid is better suitable for multilooking interferogram coordinates
        azis = np.arange(coarsen[0]//2 + 0.5, azi_max+coarsen[0] + 0.5, coarsen[0], dtype=np.float64)
        rngs = np.arange(coarsen[1]//2 + 0.5, rng_max+coarsen[1] + 0.5, coarsen[1], dtype=np.float64)

        # DEM extent in radar coordinates
        extent_ra = self.get_extent_ra()
        minx, miny, maxx, maxy = np.round(extent_ra.bounds)
        #print ('minx, miny, maxx, maxy', minx, miny, maxx, maxy)
        azis = azis[np.where((azis>=miny)&(azis<=maxy))]
        rngs = rngs[np.where((rngs>=minx)&(rngs<=maxx))]
    
        return (azis, rngs)

    def get_trans(self):
        """
        Retrieve the transform data.

        This function opens a NetCDF dataset, which contains data mapping from radar
        coordinates to geographical coordinates (from azimuth-range to latitude-longitude domain).

        Parameters
        ----------
        Returns
        -------
        xarray.Dataset or list of xarray.Dataset
            An xarray dataset(s) with the transform data.

        Examples
        --------
        Get the inverse transform data:
        get_trans()
        """
        return self.open_cube('trans')

    def compute_trans(self, coarsen, dem='auto', interactive=False):
        """
        Retrieve or calculate the transform data. This transform data is then saved as
        a NetCDF file for future use.

        This function generates data mapping from geographical coordinates to radar coordinates (azimuth-range domain).
        The function uses a Digital Elevation Model (DEM) to derive the geographical coordinates, and then uses the
        `SAT_llt2rat` function to map these to radar coordinates.

        Parameters
        ----------
        coarsen : int or (int, int)
            The decimation factor in the azimuth and range direction.

        Returns
        -------
        None

        Examples
        --------
        Calculate and get the transform data:
        >>> Stack.compute_trans_dat(1)
        """
        import dask
        import xarray as xr
        import numpy as np
        import os
        from tqdm.auto import tqdm
        import joblib
        import warnings
        warnings.filterwarnings('ignore')

        # range, azimuth, elevation(ref to radius in PRM), lon, lat [ASCII default] 
        #llt2rat_map = {0: 'rng', 1: 'azi', 2: 'ele', 3: 'll', 4: 'lt'}
        # use only 3 values from 5 available ignoring redundant lat, lon
        llt2rat_map = {0: 'rng', 1: 'azi', 2: 'ele'}

        coarsen = self.get_coarsen(coarsen)

        prm = self.PRM()
        def SAT_llt2rat(lats, lons, zs):
            # for binary=True values outside of the scene missed and the array is not complete
            # 4th and 5th coordinates are the same as input lat, lon
            return prm.SAT_llt2rat(np.column_stack([lons, lats, zs]),
                                        precise=1, binary=False)\
                           .astype(np.float32).reshape(zs.size, 5)[...,:3]

        # exclude latitude and longitude columns as redundant
        def trans_block(lats, lons, amin=-np.inf, amax=np.inf, rmin=-np.inf, rmax=np.inf, filename=None):
            # disable "distributed.utils_perf - WARNING - full garbage collections ..."
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()
            import warnings
            warnings.filterwarnings('ignore')

            dlat = dem.lat.diff('lat')[0]
            dlon = dem.lon.diff('lon')[0]
            topo = dem.sel(lat=slice(lats[0]-dlat, lats[-1]+dlat), lon=slice(lons[0]-dlon, lons[-1]+dlon))\
                      .compute(n_workers=1)
            #print ('topo.shape', topo.shape, 'lats.size', lats.size, 'lons', lons.size)
            if np.isfinite(amin):
                # check if the topo block is empty or not
                lts = topo.lat.values
                lls = topo.lon.values
                border_lts = np.concatenate([lts, lts, np.repeat(lts[0], lls.size), np.repeat(lts[-1], lls.size)])
                border_lls = np.concatenate([np.repeat(lls[0], lts.size), np.repeat(lls[-1], lts.size), lls, lls])
                border_zs  = np.concatenate([topo.values[:,0], topo.values[:,-1], topo.values[0,:], topo.values[-1,:]])
                rae = SAT_llt2rat(border_lts, border_lls, border_zs)
                del lts, lls, border_lts, border_lls, border_zs
                # this mask does not work for a single chunk
                #mask = (rae[:,0]>=rmin) & (rae[:,0]<=rmax) & (rae[:,1]>=amin) & (rae[:,1]<=amax)
                #del rae
                #valid_pixels = mask[mask].size > 0
                #del mask
                invalid_mask = ((rae[:,0]<rmin) | (rmax<rae[:,0])) & ((rae[:,1]<amin) | (amax<rae[:,1]))
                del rae
                valid_pixels = invalid_mask[~invalid_mask].size > 0
                del invalid_mask
            else:
                # continue the processing without empty block check
                valid_pixels = True

            if valid_pixels:
                grid = topo.interp({topo.dims[0]: lats, topo.dims[1]: lons})
                del topo
                # compute 3D radar coordinates for all the geographical 3D points
                lls, lts = np.meshgrid(lons.astype(np.float32), lats.astype(np.float32))
                rae = SAT_llt2rat(lts.ravel(), lls.ravel(), grid.values.ravel())
                del lls, lts, grid
                # mask invalid values for better compression
                mask = (rae[...,0]>=rmin) & (rae[...,0]<=rmax) & (rae[...,1]>=amin) & (rae[...,1]<=amax)
                rae[~mask] = np.nan
                del mask
                rae = rae.reshape(lats.size, lons.size, -1).transpose(2,0,1)
            else:
                rae = np.nan * np.zeros((3, lats.size, lons.size), np.float32)

            if filename is None:
                return rae
            # transform to separate variables, round for better compression
            trans = xr.Dataset({val: xr.DataArray(rae[key],
                            coords={'lat': lats,'lon': lons}) for (key, val) in llt2rat_map.items()})
            encoding = {vn: self._compression(trans[vn].shape, chunksize=self.netcdf_chunksize) for vn in trans.data_vars}
            if os.path.exists(filename):
                os.remove(filename)
            trans.to_netcdf(filename, encoding=encoding, engine=self.netcdf_engine)
            del trans

        if isinstance(dem, str) and dem == 'auto':
            # do not use coordinate names lat,lon because the output grid saved as (lon,lon) in this case...
            dem = self.get_dem()
        #dem = dem.rename({'lat': 'yy', 'lon': 'xx'})

        # check DEM corners
        dem_corners = dem[::dem.lat.size-1, ::dem.lon.size-1].compute()
        rngs, azis, _ = trans_block(dem_corners.lat.values, dem_corners.lon.values)
        azi_size = abs(np.diff(azis, axis=0).mean())
        rng_size = abs(np.diff(rngs, axis=1).mean())
        del rngs, azis, _
        #print ('azi_size', azi_size)
        #print ('rng_size', rng_size)
        azi_steps = int(np.round(azi_size / coarsen[0]))
        rng_steps = int(np.round(rng_size / coarsen[1]))
        #print ('azi_steps', azi_steps, 'rng_steps',rng_steps)

        # select radar coordinates extent
        azis, rngs = self.define_trans_grid(coarsen)
        azi_max = np.max(azis)
        rng_max = np.max(rngs)
        borders = {'amin': -coarsen[0], 'amax': azi_max + coarsen[0],
                   'rmin': -coarsen[1], 'rmax': rng_max + coarsen[1]}
        #print ('borders', borders)

        # process the area
        lats = np.linspace(dem.lat[0], dem.lat[-1], azi_steps)
        lons = np.linspace(dem.lon[0], dem.lon[-1], rng_steps)
        #print ('lats', lats, 'lons', lons)
        #print ('lats.size', lats.size, 'lons.size', lons.size)

        # split to equal chunks and rest
        lats_blocks = np.array_split(lats, np.arange(0,lats.size, self.chunksize)[1:])
        lons_blocks = np.array_split(lons, np.arange(0,lons.size, self.chunksize)[1:])
        #print ('lats_blocks.size', len(lats_blocks), 'lons_blocks.size', len(lons_blocks))
        #print ('lats_blocks[0]', lats_blocks[0])

        # helper function
        chunks = len(lats_blocks), len(lons_blocks)
        digits = len(str(chunks[0]*chunks[1]))
        def fullname(index):
            return os.path.join(self.basedir, f'trans_{index:0{digits}d}.grd')
        # Process data in chunks and save each chunk into a separate NetCDF file.
        with self.tqdm_joblib(tqdm(desc='Radar Transform Computing', total=chunks[0]*chunks[1])) as progress_bar:
            joblib.Parallel(n_jobs=-1)(joblib.delayed(trans_block)(
                lats_blocks[lat], lons_blocks[lon], **borders, filename=fullname(index)
            ) for index, (lat, lon) in enumerate(np.ndindex(chunks[0], chunks[1])))

        filenames = [fullname(index[0]) for index in enumerate(np.ndindex(chunks[0], chunks[1]))]
        # re-save all the chunk NetCDF files as single NetCDF file
        trans = xr.open_mfdataset(
            np.asarray(filenames).reshape((chunks[0], chunks[1])).tolist(),
            engine=self.netcdf_engine,
            chunks=self.chunksize,
            parallel=True,
            concat_dim=['lat','lon'],
            combine='nested'
        )
        # fix geographic coordinates
        #print ('lats', np.diff(lats)[:10])
        #print ('trans.lat', np.diff(trans.lat)[10])
        # add target radar coordinate grid for the user defined spacing (coarsen)
        azis, rngs = self.define_trans_grid(coarsen)
        trans['y'] = azis
        trans['x'] = rngs

        if interactive:
            return trans
        else:
            # use safe=False attribute to save y,x coordinates as is
            self.save_cube(trans, 'trans', 'Radar Transform Saving')
            del lats_blocks, lons_blocks, trans
            # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
            import gc; gc.collect()
        # cleanup
        for filename in filenames:
            os.remove(filename)
        del filenames
