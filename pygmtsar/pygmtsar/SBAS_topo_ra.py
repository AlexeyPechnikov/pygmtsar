# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_trans import SBAS_trans
from .PRM import PRM
from .tqdm_dask import tqdm_dask

class SBAS_topo_ra(SBAS_trans):

    def topo_ra(self, subswath=None, coarsen=(2,2), chunksize=None, interactive=False):
        """
        Compute the topography in radar coordinates (topo_ra).

        Parameters
        ----------
        subswath : int or None, optional
            The subswath number to compute the topographic radar coordinates for. If None, the computation
            will be performed for all subswaths. Default is None.
        idec : int, optional
            The decimation factor in the range direction. Default is 2.
        jdec : int, optional
            The decimation factor in the azimuth direction. Default is 2.
        n_jobs : int, optional
            The number of parallel jobs to run. Default is -1, which uses all available CPU cores.
        interactive : bool, optional
            If True, the computation will be performed interactively and the result will be returned as a delayed object.
            Default is False.

        Notes
        -----
        This method computes the topography in radar coordinates (topo_ra) for the specified subswath(s). It uses the trans.dat
        and trans_blocks_extents data files to build the necessary index tree and perform the coordinate transformation.
        The computed topo_ra grids are saved in NetCDF files and can be accessed using the 'get_topo_ra' method.
        """
        from scipy.spatial import cKDTree
        import dask
        import xarray as xr
        import numpy as np
        import os

#         # extract and process a single trans_dat subset
#         @dask.delayed
#         def topo_ra_block_prepare(azis, rngs, chunksize=None):
#             dazi = np.diff(azis)[0]
#             drng = np.diff(rngs)[0]
#             azis_min = azis.min() - dazi
#             azis_max = azis.max() + dazi
#             rngs_min = rngs.min() - drng
#             rngs_max = rngs.max() + drng
#             #print ('azis_min', azis_min, 'azis_max', azis_max)
# 
#             lats = trans_dat.lat[((trans_dat.azi_min<=azis_max)&(trans_dat.azi_max>=azis_min))]
#             lons = trans_dat.lon[((trans_dat.rng_min<=rngs_max)&(trans_dat.rng_max>=rngs_min))]
#             #print ('lats.shape', lats.shape, 'lons.shape', lons.shape)
# 
#             # extract and materialize required subset
#             trans_subset = trans_dat.sel(lat=lats, lon=lons)
#             block_ele = trans_subset.ele.compute(n_workers=1).data.ravel()
#             block_azi = trans_subset.azi.compute(n_workers=1).data.ravel()
#             block_rng = trans_subset.rng.compute(n_workers=1).data.ravel()
#             mask = (block_azi>=azis_min)&(block_azi<=azis_max)&(block_rng>=rngs_min)&(block_rng<=rngs_max)
#             block_ele_masked = block_ele[mask]
#             block_azi_masked = block_azi[mask]
#             block_rng_masked = block_rng[mask]
# 
#             del lats, lons, trans_subset, mask
#             return (block_azi, block_rng, block_ele)
        
        # extract and process multiple chunked trans_dat subsets
        # it can be some times slow and requires much less memory
        @dask.delayed
        def topo_ra_block_prepare(azis, rngs, chunksize):
            dazi = np.diff(azis)[0]
            drng = np.diff(rngs)[0]
            azis_min = azis.min() - dazi
            azis_max = azis.max() + dazi
            rngs_min = rngs.min() - drng
            rngs_max = rngs.max() + drng
            #print ('azis_min', azis_min, 'azis_max', azis_max)

            lats = trans_dat.lat[((trans_dat.azi_min<=azis_max)&(trans_dat.azi_max>=azis_min))]
            lons = trans_dat.lon[((trans_dat.rng_min<=rngs_max)&(trans_dat.rng_max>=rngs_min))]
            #print ('lats.shape', lats.shape, 'lons.shape', lons.shape)

            # split to equal chunks and rest
            blocks = int(np.ceil(lats.size*lons.size / chunksize**2))
            lats_blocks = np.array_split(lats, np.arange(0, lats.size, chunksize**2 // lons.size)[1:])
            # process chunks
            block_azis = []
            block_rngs = []
            block_eles = []
            for lats_block in lats_blocks:
                # extract and materialize required subset
                trans_subset = trans_dat.sel(lat=lats_block, lon=lons)
                block_ele = trans_subset.ele.compute(n_workers=1).data.ravel()
                block_azi = trans_subset.azi.compute(n_workers=1).data.ravel()
                block_rng = trans_subset.rng.compute(n_workers=1).data.ravel()
                mask = (block_azi>=azis_min)&(block_azi<=azis_max)&(block_rng>=rngs_min)&(block_rng<=rngs_max)
                block_eles.append(block_ele[mask])
                block_azis.append(block_azi[mask])
                block_rngs.append(block_rng[mask])
                # cleanup
                del trans_subset, block_ele, block_azi, block_rng, mask
            # merge extracted results
            block_azi = np.concatenate(block_azis)
            block_rng = np.concatenate(block_rngs)
            block_ele = np.concatenate(block_eles)

            del block_azis, block_rngs, block_eles
            return (block_azi, block_rng, block_ele)

#         # cKDTree interpolations allows to get the distances to nearest pixels
#         @dask.delayed
#         def topo_ra_block(data, azis, rngs):
#             from scipy.spatial import cKDTree
# 
#             block_azi, block_rng, block_ele = data
# 
#             # interpolate topo_ra on trans_dat
#             grid_azi, grid_rng = np.meshgrid(azis, rngs)
#             tree = cKDTree(np.column_stack([block_azi, block_rng]), compact_nodes=False, balanced_tree=False)
#             d, inds = tree.query(np.column_stack([grid_azi.ravel(), grid_rng.ravel()]), k = 1, workers=1)
#             grid = block_ele[inds]
#             #print ('distance range', d.min().round(2), d.max().round(2))
# 
#             del block_ele, block_azi, block_rng, tree
#             return grid.reshape((rngs.size, azis.size)).T

        # griddata interpolation is easy and provides multiple methods
        @dask.delayed
        def topo_ra_block(data, azis, rngs):
            from scipy.interpolate import griddata

            block_azi, block_rng, block_ele = data
            points = np.column_stack([block_azi, block_rng])

            # interpolate topo_ra on trans_dat
            grid = griddata(points, block_ele, (azis[None, :], rngs[:, None]), method='nearest')
            del block_ele, block_azi, block_rng, points
            return grid.reshape((rngs.size, azis.size)).T

        if chunksize is None:
            chunksize = self.chunksize

        # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        trans_dat = self.get_trans_dat(subswath)
        # materialize indices
        trans_dat['azi_min'] = trans_dat.azi_min.compute()
        trans_dat['azi_max'] = trans_dat.azi_max.compute()
        trans_dat['rng_min'] = trans_dat.rng_min.compute()
        trans_dat['rng_max'] = trans_dat.rng_max.compute()

        # define topo_ra grid
        rng_max, yvalid, num_patch = self.PRM(subswath).get('num_rng_bins', 'num_valid_az', 'num_patches')
        azi_max = yvalid * num_patch
        #print ('DEBUG: rng_max', rng_max, 'azi_max', azi_max)
        # use center pixel GMT registration mode
        azis = np.arange(1, azi_max+1, coarsen[0], dtype=np.int32)
        rngs = np.arange(1, rng_max+1, coarsen[1], dtype=np.int32)

        # build topo_ra grid by chunks

        # split to equal chunks and rest
        azis_blocks = np.array_split(azis, np.arange(0, azis.size, chunksize)[1:])
        rngs_blocks = np.array_split(rngs, np.arange(0, rngs.size, chunksize)[1:])
        #print ('azis_blocks.size', len(azis_blocks), 'rngs_blocks.size', len(rngs_blocks))
        #print ('lats_blocks[0]', lats_blocks[0])

    #     # single block test
    #     data = topo_ra_block_prepare(azis_blocks[-1], rngs_blocks[-1])
    #     block = topo_ra_block(data, azis_blocks[-1], rngs_blocks[-1])
    #     print (block.compute())

    #     # single block test
    #     #print ('azis_blocks[-1]', azis_blocks[-1])
    #     #print ('rngs_blocks[5]', rngs_blocks[5])
    #     data = topo_ra_block_prepare(azis_blocks[-1], rngs_blocks[5])
    #     block = topo_ra_block(data, azis_blocks[-1], rngs_blocks[5])
    #     return np.flipud(block.compute())

        blocks_total = []
        for rngs_block in rngs_blocks:
            blocks = []
            for azis_block in azis_blocks:
                # extract multiple outputs
                #blockset = dask.delayed(trans_block)(...)
                #block = dask.array.from_delayed(blockset[0], shape=(...), dtype=np.float32)
                data = topo_ra_block_prepare(azis_block, rngs_block, chunksize)
                block = dask.array.from_delayed(topo_ra_block(data, azis_block, rngs_block),
                                                shape=(azis_block.size, rngs_block.size), dtype=np.float32)
                blocks.append(block.transpose(1,0))
            blocks_total.append(blocks)
        topo_ra = dask.array.block(blocks_total).transpose(1,0)
        topo_ra = xr.DataArray(topo_ra, coords={'y': azis, 'x': rngs}).rename('topo_ra')

        if interactive:
            # do not flip vertically because it's returned as is without SBAS.get_topo_ra() function
            return topo_ra

        # save to NetCDF file
        filename = self.get_filenames(subswath, None, 'topo_ra')
        if os.path.exists(filename):
            os.remove(filename)
        # flip vertically for GMTSAR compatibility reasons
        topo_ra = xr.DataArray(dask.array.flipud(topo_ra), coords=topo_ra.coords, name=topo_ra.name)
        # rename to save lazy NetCDF preventing broken coordinates (y,y) 
        topo_ra = topo_ra.rename({'y': 'a', 'x': 'r'})
        handler = topo_ra.to_netcdf(filename,
                                    encoding={'topo_ra': self.compression(chunksize=chunksize)},
                                    engine=self.engine,
                                    compute=False)
        return handler

    def topo_ra_parallel(self, interactive=False, **kwargs):
        """
        Build topography in radar coordinates from WGS84 DEM using parallel computation.

        Parameters
        ----------
        interactive : bool, optional
            If True, the computation will be performed interactively and the results will be returned as delayed objects.
            If False, the progress will be displayed using tqdm_dask. Default is False.

        Returns
        -------
        handler or list of handlers
            The handler(s) of the delayed computation if 'interactive' is True. Otherwise, None.

        Examples
        --------
        sbas.topo_ra_parallel()

        Notes
        -----
        This method performs the parallel computation of topography in the radar coordinates (topo_ra) for all subswaths
        using Dask. It calls the 'topo_ra' method for each subswath in parallel. If 'interactive' is True, the delayed
        computation handlers will be returned. Otherwise, the progress will be displayed using tqdm_dask.
        """
        import dask

        # auto generate the trans.dat file
        self.trans_dat_parallel(**kwargs)

        # process all the subswaths
        subswaths = self.get_subswaths()
        delayeds = []
        for subswath in subswaths:
            delayed = self.topo_ra(subswath=subswath, interactive=interactive, **kwargs)
            if not interactive:
                tqdm_dask(dask.persist(delayed), desc=f'Radar Topography Computing sw{subswath}')
            else:
                delayeds.append(delayed)

        if interactive:
            return delayeds[0] if len(delayeds)==1 else delayeds

    def get_topo_ra(self):
        """
        Get the radar topography grid.

        Returns
        -------
        xarray.DataArray or list of xarray.DataArray
            The 'topo_ra' grid data as a single xarray.DataArray if only one grid is found,
            or a list of xarray.DataArray if multiple grids are found.

        Examples
        --------
        Get DEM for all the processed subswaths:
        topo_ra = sbas.get_topo_ra()

        Notes
        -----
        This method opens the 'topo_ra' grids using the `open_grids` method and applies the `func` function
        to each grid to flip it vertically for compatibility reasons with GMTSAR. The resulting grids are returned.
        """
        import xarray as xr
        import dask.array

        def func(topo):
            # flip vertically for GMTSAR compatibility reasons
            topo = xr.DataArray(dask.array.flipud(topo), coords=topo.coords, attrs=topo.attrs, name=topo.name)
            # fix renamed dimensions required to save lazy NetCDF properly
            if 'a' in topo.dims and 'r' in topo.dims:
                return topo.rename({'a': 'y', 'r': 'x'})
            return topo

        topos = self.open_grids(None, 'topo_ra', func=func)

        return topos[0] if len(topos)==1 else topos
