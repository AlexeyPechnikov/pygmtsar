# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_stack import SBAS_stack
from .tqdm_dask import tqdm_dask

class SBAS_trans(SBAS_stack):
    
    def define_trans_grid(self, subswath, coarsen):
        import numpy as np
        # select radar coordinates extent
        rng_max, yvalid, num_patch = self.PRM(subswath).get('num_rng_bins', 'num_valid_az', 'num_patches')
        azi_max = yvalid * num_patch
        #print ('azi_max', azi_max, 'rng_max', rng_max)
        # this grid covers the full interferogram area
        azis = np.arange(0, azi_max+coarsen[0], coarsen[0], dtype=np.int32)
        rngs = np.arange(0, rng_max+coarsen[1], coarsen[1], dtype=np.int32)
        # this grid is better suitable for multilooking interferogram coordinates
        #azis = np.arange(coarsen[0]//2, azi_max+coarsen[0], coarsen[0], dtype=np.int32)
        #rngs = np.arange(coarsen[1]//2, rng_max+coarsen[1], coarsen[1], dtype=np.int32)
        return (azis, rngs)
        
    def get_trans_dat(self, subswath=None, chunksize=None):
        """
        Retrieve the transform data for a specific or all subswaths.

        This function opens a NetCDF dataset, which contains data mapping from geographical
        coordinates to radar coordinates (from latitude-longitude domain to azimuth-range).

        Parameters
        ----------
        subswath : int, optional
            Subswath number to retrieve. If not specified, the function will retrieve the transform
            data for all available subswaths.

        Returns
        -------
        xarray.Dataset
            An xarray dataset with the transform data.

        Examples
        --------
        Get the transform data for a specific subswath:
        get_trans_dat(1)

        Get the transform data for all available subswaths:
        get_trans_dat()
        """
        import xarray as xr
    
        if subswath is None:
            subswaths = self.get_subswaths()
        else:
            subswaths = [subswath]

        if chunksize is None:
            chunksize = self.chunksize

        transs = []
        for subswath in subswaths:
            filename = self.get_filenames(subswath, None, 'trans')
            trans = xr.open_dataset(filename, engine=self.engine, chunks=chunksize)
            if 'yy' in trans and 'xx' in trans:
                transs.append(trans.rename({'yy': 'lat', 'xx': 'lon'}))
            else:
                transs.append(trans)
        return transs[0] if len(transs)==1 else transs

    def trans_dat(self, subswath=None, coarsen=2, chunksize=None, interactive=False):
        """
        Retrieve or calculate the transform data for a specific or all subswaths. This transform data is then saved as
        a NetCDF file for future use.

        This function generates data mapping from geographical coordinates to radar coordinates (azimuth-range domain).
        The function uses a Digital Elevation Model (DEM) to derive the geographical coordinates, and then uses the
        `SAT_llt2rat` function to map these to radar coordinates.

        Parameters
        ----------
        subswath : int, optional
            Subswath number to retrieve. If not specified, the function will retrieve the transform
            data for all available subswaths.
        coarsen(jdec, idec) : (int, int) , optional
            The decimation factor in the azimuth and range direction. Default is 2.
        interactive : bool, optional
            If True, the function returns the transform data without saving it. If False, the function
            saves the transform data as a NetCDF file. Default is False.

        Returns
        -------
        xarray.Dataset or dask.delayed.Delayed
            If interactive is True, it returns an xarray dataset with the transform data.
            If interactive is False, it returns a dask Delayed object representing the computation of writing
            the transform data to a NetCDF file.

        Examples
        --------
        Calculate and get the transform data for a specific subswath:
        >>> trans_dat(1)
        <dask.delayed.Delayed at 0x7f8d13a69a90>

        Calculate and get the transform data for all available subswaths:
        >>> trans_dat()
        <dask.delayed.Delayed at 0x7f8d13a69b70>

        Calculate and get the transform data without saving it:
        >>> trans_dat(interactive=True)
        """
        import dask
        import xarray as xr
        import numpy as np
        import os
        import sys

        # range, azimuth, elevation(ref to radius in PRM), lon, lat [ASCII default] 
        #llt2rat_map = {0: 'rng', 1: 'azi', 2: 'ele', 3: 'll', 4: 'lt'}
        # use only 3 values from 5 available ignoring redudant lat, lon
        llt2rat_map = {0: 'rng', 1: 'azi', 2: 'ele'}

        if chunksize is None:
            chunksize = self.chunksize

        # expand simplified definition
        if np.issubdtype(type(coarsen), np.integer):
            coarsen = (coarsen, coarsen)

        # build trans.dat
        def SAT_llt2rat(z, lat, lon, subswath, amin=0, amax=None, rmin=0, rmax=None):
            if amax is not None and rmax is not None:
                # check border coordinates to detect if the block is completely outside of radar area
                lats = np.concatenate([lat, lat, np.repeat(lat[0], lon.size), np.repeat(lat[-1], lon.size)])
                lons = np.concatenate([np.repeat(lon[0], lat.size), np.repeat(lon[-1], lat.size), lon, lon])
                zs = np.concatenate([z[:,0], z[:,-1], z[0,:], z[-1,:]])
                coords_ll = np.column_stack([lons, lats, zs])
                # raell
                coords_ra = self.PRM(subswath).SAT_llt2rat(coords_ll, precise=1, binary=False)\
                    .astype(np.float32).reshape(zs.size, 5)
                #print (coords_ra.shape)
                # check validity
                mask = (coords_ra[:,0]>=rmin) & (coords_ra[:,0]<=rmax) & (coords_ra[:,1]>=amin) & (coords_ra[:,1]<=amax)
                #print ('mask[mask].size', mask[mask].size)
                if mask[mask].size == 0:
                    return np.nan * np.zeros((z.shape[0], z.shape[1], 5), np.float32)

            # compute 3D radar coordinates for all the geographical 3D points
            lons, lats = np.meshgrid(lon.astype(np.float32), lat.astype(np.float32))
            coords_ll = np.column_stack([lons.ravel(), lats.ravel(), z.ravel()])
            # for binary=True values outside of the scene missed and the array is not complete
            # 4th and 5th coordinates are the same as input lat, lon
            coords_ra = self.PRM(subswath).SAT_llt2rat(coords_ll, precise=1, binary=False)\
                .astype(np.float32).reshape(z.shape[0], z.shape[1], 5)[...,:3]
            if amax is not None and rmax is not None:
                # mask values outside of radar area
                mask = (coords_ra[...,0]>=rmin) & (coords_ra[...,0]<=rmax) & (coords_ra[...,1]>=amin) & (coords_ra[...,1]<=amax)
                coords_ra[~mask] = np.nan
            return coords_ra

        # exclude latitude and longitude columns as redudant
        def trans_block(lats, lons, subswath, amin=0, amax=None, rmin=0, rmax=None):
            dlat = dem.yy.diff('yy')[0]
            dlon = dem.xx.diff('xx')[0]
            topo = dem.sel(yy=slice(lats[0]-dlat, lats[-1]+dlat), xx=slice(lons[0]-dlon, lons[-1]+dlon))
            #print ('topo.shape', topo.shape, 'lats.size', lats.size, 'lons', lons.size)
            grid = topo.interp({topo.dims[0]: lats, topo.dims[1]: lons})
            #print ('grid.shape', grid.shape) 
            rae = SAT_llt2rat(grid.values, lats, lons, subswath, amin, amax, rmin, rmax)
            # define radar coordinate extent for the block
            # this code produces a lot of warnings "RuntimeWarning: All-NaN slice encountered"
            #rmin, rmax = np.nanmin(rae[...,0]), np.nanmax(rae[...,0])
            #amin, amax = np.nanmin(rae[...,1]), np.nanmax(rae[...,1])
            # this code spends additional time for the checks to exclude warnings
            if not np.all(np.isnan(rae[...,0])):
                rmin, rmax = np.nanmin(rae[...,0]), np.nanmax(rae[...,0])
            else:
                rmin = rmax = np.nan
            if not np.all(np.isnan(rae[...,1])):
                amin, amax = np.nanmin(rae[...,1]), np.nanmax(rae[...,1])
            else:
                amin = amax = np.nan
            extent = np.asarray([amin, amax, rmin, rmax], dtype=np.float32)
            # cleanup to fix unmanaged memory increasing due to external object use
            del topo, grid
            return (rae, extent)

        # do not use coordinate names lat,lon because the output grid saved as (lon,lon) in this case...
        dem = self.get_dem(geoloc=True).rename({'lat': 'yy', 'lon': 'xx'})

        # check DEM corners
        dem_corners = dem[::dem.yy.size-1, ::dem.xx.size-1]
        rngs, azis, _ = trans_block(dem_corners.yy, dem_corners.xx, subswath)[0].transpose(2,0,1)
        azi_size = abs(np.diff(azis, axis=0).mean())
        rng_size = abs(np.diff(rngs, axis=1).mean())
        #print ('azi_size', azi_size)
        #print ('rng_size', rng_size)
        azi_steps = int(np.round(azi_size / coarsen[0]))
        rng_steps = int(np.round(rng_size / coarsen[1]))
        #print ('azi_steps', azi_steps, 'rng_steps',rng_steps)

        # select radar coordinates extent
        azis, rngs = self.define_trans_grid(subswath, coarsen)
        azi_max = np.max(azis)
        rng_max = np.max(rngs)
        # allow 2 points around for linear and cubic interpolations
        borders = {'amin': -2*azi_size/azi_steps, 'amax': azi_max+2*azi_size/azi_steps,
                 'rmin': -2*rng_size/rng_steps, 'rmax': rng_max+2*rng_size/rng_steps}
        #print ('borders', borders)

        ################################################################################
        # process the area
        ################################################################################
        lats = np.linspace(dem.yy[0], dem.yy[-1], azi_steps)
        lons = np.linspace(dem.xx[0], dem.xx[-1], rng_steps)
        #print ('lats', lats, 'lons', lons)
        #print ('lats.size', lats.size, 'lons.size', lons.size)

        # split to equal chunks and rest
        lats_blocks = np.array_split(lats, np.arange(0,lats.size, chunksize)[1:])
        lons_blocks = np.array_split(lons, np.arange(0,lons.size, chunksize)[1:])
        #print ('lats_blocks.size', len(lats_blocks), 'lons_blocks.size', len(lons_blocks))
        #print ('lats_blocks[0]', lats_blocks[0])

        blocks_total = []
        extents_total = []
        for lons_block in lons_blocks:
            blocks = []
            extents = []
            for lats_block in lats_blocks:
                # extract multiple outputs
                blockset = dask.delayed(trans_block)(lats_block, lons_block, subswath, **borders)
                block = dask.array.from_delayed(blockset[0], shape=(lats_block.size, lons_block.size, 3), dtype=np.float32)
                extent = dask.array.from_delayed(blockset[1], shape=(4,), dtype=np.float32)
                #block = dask.array.from_delayed(dask.delayed(trans_block)(lats_block, lons_block, subswath, **extent),
                #                                shape=(lats_block.size, lons_block.size, 3), dtype=np.float32)
                blocks.append(block.transpose(2,1,0))
                extents.append(extent[:, None, None])
                del blockset, block, extent
            blocks_total.append(blocks)
            extents_total.append(extents)
            del blocks, extents
        rae = dask.array.block(blocks_total).transpose(2,1,0)
        extent = dask.array.block(extents_total).transpose(2,1,0)
        del blocks_total, extents_total

        # transform to separate variables
        key_boundaries = {val: xr.DataArray(extent[...,key], dims=['block_azi', 'block_rng'])
                          for (key, val) in enumerate(['amin', 'amax', 'rmin', 'rmax'])}
        keys_vars = {val: xr.DataArray(rae[...,key], coords={'lat': lats,'lon': lons})
                          for (key, val) in llt2rat_map.items()}
        trans = xr.Dataset({**key_boundaries, **keys_vars})

        # add corresponding radar grid coordinates
        trans['y'] = azis
        trans['x'] = rngs

        if interactive:
            return trans

        # save to NetCDF file
        filename = self.get_filenames(subswath, None, 'trans')
        if os.path.exists(filename):
            os.remove(filename)
        encoding = {val: self.compression(trans[val].shape, chunksize=chunksize) for (key, val) in llt2rat_map.items()}
        handler = trans.to_netcdf(filename,
                                        encoding=encoding,
                                        engine=self.engine,
                                        compute=False)
        return handler

    def trans_dat_parallel(self, interactive=False, **kwargs):
        """
        Retrieve or calculate the transform data for all subswaths in parallel. This function processes each subswath
        concurrently using Dask.

        Parameters
        ----------
        interactive : bool, optional
            If True, the function returns a list of dask.delayed.Delayed objects representing the computation of
            transform data for each subswath. If False, the function processes the transform data for each subswath
            concurrently using Dask. Default is False.

        Returns
        -------
        list or dask.delayed.Delayed
            If interactive is True, it returns a list of dask.delayed.Delayed objects representing the computation
            of transform data for each subswath.
            If interactive is False, it returns a dask.delayed.Delayed object representing the computation of
            transform data for all subswaths.

        Examples
        --------
        Calculate and get the transform data for all subswaths in parallel:

        >>> trans_dat_parallel()
        <dask.delayed.Delayed at 0x7f8d13a69a90>

        Calculate and get the transform data for all subswaths in parallel without saving it:

        >>> trans_dat_parallel(interactive=True)
        [<dask.delayed.Delayed at 0x7f8d13a69a90>, <dask.delayed.Delayed at 0x7f8d13a69b70>]
        """
        import dask

        # process all the subswaths
        subswaths = self.get_subswaths()
        delayeds = []
        for subswath in subswaths:
            delayed = self.trans_dat(subswath=subswath, interactive=interactive, **kwargs)
            if not interactive:
                tqdm_dask(dask.persist(delayed), desc=f'Radar Transform Computing sw{subswath}')
            else:
                delayeds.append(delayed)

        if interactive:
            return delayeds[0] if len(delayeds)==1 else delayeds

        # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
        import gc; gc.collect()
