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
        azis = np.arange(coarsen[0]//2, azi_max+coarsen[0], coarsen[0], dtype=np.int32)
        rngs = np.arange(coarsen[1]//2, rng_max+coarsen[1], coarsen[1], dtype=np.int32)
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
        return self.open_grid('trans')

    def trans(self, coarsen, buffer_degrees=None, interactive=False):
        """
        Retrieve or calculate the transform data. This transform data is then saved as
        a NetCDF file for future use.

        This function generates data mapping from geographical coordinates to radar coordinates (azimuth-range domain).
        The function uses a Digital Elevation Model (DEM) to derive the geographical coordinates, and then uses the
        `SAT_llt2rat` function to map these to radar coordinates.

        Parameters
        ----------
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
        Calculate and get the transform data:
        >>> trans_dat()
        <dask.delayed.Delayed at 0x7f8d13a69a90>

        Calculate and get the transform data without saving it:
        >>> trans_dat(interactive=True)
        """
        import dask
        import xarray as xr
        import numpy as np
        #import os
        #import sys

        # range, azimuth, elevation(ref to radius in PRM), lon, lat [ASCII default] 
        #llt2rat_map = {0: 'rng', 1: 'azi', 2: 'ele', 3: 'll', 4: 'lt'}
        # use only 3 values from 5 available ignoring redudant lat, lon
        llt2rat_map = {0: 'rng', 1: 'azi', 2: 'ele'}

        # expand simplified definition
        if not isinstance(coarsen, (list,tuple, np.ndarray)):
            coarsen = (coarsen, coarsen)

        # build trans.dat
        def SAT_llt2rat(z, lat, lon, amin=0, amax=None, rmin=0, rmax=None):
            if amax is not None and rmax is not None:
                # check border coordinates to detect if the block is completely outside of radar area
                lats = np.concatenate([lat, lat, np.repeat(lat[0], lon.size), np.repeat(lat[-1], lon.size)])
                lons = np.concatenate([np.repeat(lon[0], lat.size), np.repeat(lon[-1], lat.size), lon, lon])
                zs = np.concatenate([z[:,0], z[:,-1], z[0,:], z[-1,:]])
                coords_ll = np.column_stack([lons, lats, zs])
                # raell
                coords_ra = self.PRM().SAT_llt2rat(coords_ll, precise=1, binary=False)\
                    .astype(np.float32).reshape(zs.size, 5)
                del lons, lats, coords_ll
                #print (coords_ra.shape)
                # check validity
                mask = (coords_ra[:,0]>=rmin) & (coords_ra[:,0]<=rmax) & (coords_ra[:,1]>=amin) & (coords_ra[:,1]<=amax)
                del coords_ra
                #print ('mask[mask].size', mask[mask].size)
                valid_pixels = mask[mask].size
                del mask
                if valid_pixels == 0:
                    # no valid pixels in the block, miss the processing
                    return np.nan * np.zeros((z.shape[0], z.shape[1], 5), np.float32)
                # valid points included into the block, continue

            # compute 3D radar coordinates for all the geographical 3D points
            lons, lats = np.meshgrid(lon.astype(np.float32), lat.astype(np.float32))
            coords_ll = np.column_stack([lons.ravel(), lats.ravel(), z.ravel()])
            # for binary=True values outside of the scene missed and the array is not complete
            # 4th and 5th coordinates are the same as input lat, lon
            coords_ra = self.PRM().SAT_llt2rat(coords_ll, precise=1, binary=False).astype(np.float32)
            if coords_ra.size == 0:
                return coords_ra
            coords_ra = coords_ra.reshape(z.shape[0], z.shape[1], 5)[...,:3]
            del lons, lats, coords_ll
            if amax is not None and rmax is not None:
                # mask values outside of radar area
                mask = (coords_ra[...,0]>=rmin) & (coords_ra[...,0]<=rmax) & (coords_ra[...,1]>=amin) & (coords_ra[...,1]<=amax)
                coords_ra[~mask] = np.nan
                del mask
            return coords_ra

        # exclude latitude and longitude columns as redudant
        def trans_block(lats, lons, amin=0, amax=None, rmin=0, rmax=None):
            dlat = dem.yy.diff('yy')[0]
            dlon = dem.xx.diff('xx')[0]
            topo = dem.sel(yy=slice(lats[0]-dlat, lats[-1]+dlat), xx=slice(lons[0]-dlon, lons[-1]+dlon))
            del dlat, dlon
            #print ('topo.shape', topo.shape, 'lats.size', lats.size, 'lons', lons.size)
            grid = topo.interp({topo.dims[0]: lats, topo.dims[1]: lons})
            del topo
            #print ('grid.shape', grid.shape) 
            rae = SAT_llt2rat(grid.values, lats, lons, amin, amax, rmin, rmax)
            del grid
            # define radar coordinate extent for the block
            # this code produces a lot of warnings "RuntimeWarning: All-NaN slice encountered"
            #rmin, rmax = np.nanmin(rae[...,0]), np.nanmax(rae[...,0])
            #amin, amax = np.nanmin(rae[...,1]), np.nanmax(rae[...,1])
            # this code spends additional time for the checks to exclude warnings
            if not np.all(np.isnan(rae[...,:2])):
                extent = np.asarray([np.nanmin(rae[...,1]), np.nanmax(rae[...,1]), np.nanmin(rae[...,0]), np.nanmax(rae[...,0])], dtype=np.float32)
            else:
                extent = np.asarray([np.nan, np.nan, np.nan, np.nan], dtype=np.float32)
            return (rae, extent)

        # do not use coordinate names lat,lon because the output grid saved as (lon,lon) in this case...
        dem = self.get_dem(geoloc=True, buffer_degrees=buffer_degrees).rename({'lat': 'yy', 'lon': 'xx'})

        # check DEM corners
        dem_corners = dem[::dem.yy.size-1, ::dem.xx.size-1]
        rngs, azis, _ = trans_block(dem_corners.yy, dem_corners.xx)[0].transpose(2,0,1)
        azi_size = abs(np.diff(azis, axis=0).mean())
        rng_size = abs(np.diff(rngs, axis=1).mean())
        #print ('azi_size', azi_size)
        #print ('rng_size', rng_size)
        azi_steps = int(np.round(azi_size / coarsen[0]))
        rng_steps = int(np.round(rng_size / coarsen[1]))
        #print ('azi_steps', azi_steps, 'rng_steps',rng_steps)

        # select radar coordinates extent
        azis, rngs = self.define_trans_grid(coarsen)
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
        lats_blocks = np.array_split(lats, np.arange(0,lats.size, self.chunksize)[1:])
        lons_blocks = np.array_split(lons, np.arange(0,lons.size, self.chunksize)[1:])
        #print ('lats_blocks.size', len(lats_blocks), 'lons_blocks.size', len(lons_blocks))
        #print ('lats_blocks[0]', lats_blocks[0])

        blocks_total = []
        extents_total = []
        for lons_block in lons_blocks:
            blocks = []
            extents = []
            for lats_block in lats_blocks:
                # extract multiple outputs
                blockset = dask.delayed(trans_block)(lats_block, lons_block, **borders)
                block = dask.array.from_delayed(blockset[0], shape=(lats_block.size, lons_block.size, 3), dtype=np.float32)
                extent = dask.array.from_delayed(blockset[1], shape=(4,), dtype=np.float32)
                #block = dask.array.from_delayed(dask.delayed(trans_block)(lats_block, lons_block, **extent),
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
        del extent, rae
        trans = xr.Dataset({**key_boundaries, **keys_vars})
        del key_boundaries, keys_vars

        # add corresponding radar grid coordinates
        trans['y'] = azis
        trans['x'] = rngs
        del azis, rngs

        if interactive:
            return trans
        return self.save_grid(trans, 'trans', 'Radar Transform Computing')
