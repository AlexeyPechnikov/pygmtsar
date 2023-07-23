# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_trans import SBAS_trans
from .tqdm_dask import tqdm_dask

class SBAS_trans_inv(SBAS_trans):
    
    def get_trans_dat_inv(self, subswath=None):
        """
        Retrieve the inverse transform data for a specific or all subswaths.

        This function opens a NetCDF dataset, which contains data mapping from radar
        coordinates to geographical coordinates (from azimuth-range to latitude-longitude domain).

        Parameters
        ----------
        subswath : int, optional
            Subswath number to retrieve. If not specified, the function will retrieve the inverse
            transform data for all available subswaths.

        Returns
        -------
        xarray.Dataset
            An xarray dataset with the transform data.

        Examples
        --------
        Get the inverse transform data for a specific subswath:
        get_trans_dat_inv(1)

        Get the inverse transform data for all available subswaths:
        get_trans_dat_inv()
        """
        import xarray as xr

        subswath = self.get_subswath(subswath)
        filename = self.get_filenames(subswath, None, 'trans_inv')
        trans_inv = xr.open_dataset(filename, engine=self.engine, chunks=self.chunksize).rename({'a': 'y', 'r': 'x'})
        return trans_inv

    def trans_dat_inv(self, subswath=None, coarsen=(2,2), chunksize=None, interactive=False):
        """
        Retrieve or calculate the transform data for a specific or all subswaths. This transform data is then saved as
            a NetCDF file for future use.

            This function generates data mapping from radar coordinates to geographical coordinates.
            The function uses the direct transform data.

        Parameters
        ----------
        subswath : int or None, optional
            The subswath number to compute the topographic radar coordinates for. If None, the computation
            will be performed for all subswaths. Default is None.
        coarsen(jdec, idec) : (int, int) , optional
            The decimation factor in the azimuth and range direction. Default is 2.
        interactive : bool, optional
            If True, the computation will be performed interactively and the result will be returned as a delayed object.
            Default is False.
        """
        import dask
        import xarray as xr
        import numpy as np
        import os

        # expand simplified definition
        if np.issubdtype(type(coarsen), np.integer):
            coarsen = (coarsen, coarsen)

    #         # extract and process a single trans_dat subset
    #         @dask.delayed
    #         def trans_dat_inv_block_prepare(azis, rngs, chunksize=None):
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
        def trans_dat_inv_block_prepare(azis, rngs, chunksize):
            # required one delta around for nearest interpolation and two for linear
            dazi = np.diff(azis)[0]
            drng = np.diff(rngs)[0]
            azis_min = azis.min() - 2*dazi
            azis_max = azis.max() + 2*dazi
            rngs_min = rngs.min() - 2*drng
            rngs_max = rngs.max() + 2*drng
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
            block_lts  = []
            block_lls  = []
            block_eles = []
            for lats_block in lats_blocks:
                # extract and materialize required subset
                trans_subset = trans_dat.sel(lat=lats_block, lon=lons)
                block_lt, block_ll = xr.broadcast(trans_subset.lat, trans_subset.lon)
                block_lt  = block_lt.data.ravel()
                block_ll  = block_ll.data.ravel()
                block_ele = trans_subset.ele.compute(n_workers=1).data.ravel()
                block_azi = trans_subset.azi.compute(n_workers=1).data.ravel()
                block_rng = trans_subset.rng.compute(n_workers=1).data.ravel()
                mask = (block_azi>=azis_min)&(block_azi<=azis_max)&(block_rng>=rngs_min)&(block_rng<=rngs_max)
                block_azis.append(block_azi[mask])
                block_rngs.append(block_rng[mask])
                block_lts.append(block_lt[mask])
                block_lls.append(block_ll[mask])
                block_eles.append(block_ele[mask])
                # cleanup
                del trans_subset, block_azi, block_rng, block_lt, block_ll, block_ele, mask
            # merge extracted results
            block_azi = np.concatenate(block_azis)
            block_rng = np.concatenate(block_rngs)
            block_lt  = np.concatenate(block_lts)
            block_ll  = np.concatenate(block_lls)
            block_ele = np.concatenate(block_eles)

            del block_azis, block_rngs, block_lts, block_lls, block_eles
            return (block_azi, block_rng, block_lt, block_ll, block_ele)

    #         # cKDTree interpolations allows to get the distances to nearest pixels
    #         @dask.delayed
    #         def trans_dat_inv_block(data, azis, rngs):
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
        def trans_dat_inv_block(data, azis, rngs):
            from scipy.interpolate import griddata

            block_azi, block_rng, block_lt, block_ll, block_ele = data
            points = np.column_stack([block_azi, block_rng])

            output = []
            for block in [block_lt, block_ll, block_ele]:
                grid = griddata(points, block, (azis[None, :], rngs[:, None]), method='linear').astype(np.float32)
                output.append(grid.reshape((rngs.size, azis.size)).T)
                del grid

            del block_ele, block_azi, block_rng, block_lt, block_ll, points
            return np.asarray(output)

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
        #rng_max, yvalid, num_patch = self.PRM(subswath).get('num_rng_bins', 'num_valid_az', 'num_patches')
        #azi_max = yvalid * num_patch
        #print ('DEBUG: rng_max', rng_max, 'azi_max', azi_max)
        # produce the same grid as coarsed interferogram
        #azis = np.arange(coarsen[0]//2, azi_max+coarsen[0], coarsen[0], dtype=np.int32)
        #rngs = np.arange(coarsen[1]//2, rng_max+coarsen[1], coarsen[1], dtype=np.int32)
        azis, rngs = self.define_trans_grid(subswath, coarsen)
        
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

        block_lts_total  = []
        block_lls_total  = []
        block_eles_total = []
        for rngs_block in rngs_blocks:
            block_lts  = []
            block_lls  = []
            block_eles = []
            for azis_block in azis_blocks:
                # extract multiple outputs
                #blockset = dask.delayed(trans_block)(...)
                #block = dask.array.from_delayed(blockset[0], shape=(...), dtype=np.float32)
                data = trans_dat_inv_block_prepare(azis_block, rngs_block, chunksize)
                block_lt, block_ll, block_ele = dask.array.from_delayed(trans_dat_inv_block(data, azis_block, rngs_block),
                                                shape=(3, azis_block.size, rngs_block.size), dtype=np.float32)
                block_lts.append(block_lt.transpose(1,0))
                block_lls.append(block_ll.transpose(1,0))
                block_eles.append(block_ele.transpose(1,0))
                del data, block_lt, block_ll, block_ele
            block_lts_total.append(block_lts)
            block_lls_total.append(block_lls)
            block_eles_total.append(block_eles)
            del block_lts, block_lls, block_eles

        coords = {'y': azis, 'x': rngs}
        lt = dask.array.block(block_lts_total).transpose(1,0)
        lt = xr.DataArray(lt, coords=coords)
        ll = dask.array.block(block_lls_total).transpose(1,0)
        ll = xr.DataArray(ll, coords=coords)
        ele = dask.array.block(block_eles_total).transpose(1,0)
        ele = xr.DataArray(ele, coords=coords)
        del block_lts_total, block_lls_total, block_eles_total

        trans_inv = xr.Dataset({'lt': lt, 'll': ll, 'ele': ele})
        del lt, ll, ele
    
        if interactive:
            return trans_inv

        # save to NetCDF file
        filename = self.get_filenames(subswath, None, 'trans_inv')
        if os.path.exists(filename):
            os.remove(filename)
        # rename to save lazy NetCDF preventing broken coordinates (y,y) 
        trans_inv = trans_inv.rename({'y': 'a', 'x': 'r'})
        encoding = {key: self.compression(trans_inv[key].shape, chunksize=chunksize) for key in trans_inv.data_vars if len(trans_inv[key].dims)==2}
        handler = trans_inv.to_netcdf(filename,
                                    encoding=encoding,
                                    engine=self.engine,
                                    compute=False)
        return handler

    def trans_dat_inv_parallel(self, interactive=False, **kwargs):
        """
        Retrieve or calculate the inverse transform data for all subswaths in parallel. This function processes each subswath
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
        Calculate and get the inverse transform data for all subswaths in parallel:

        >>> trans_dat_inv_parallel()
        <dask.delayed.Delayed at 0x7f8d13a69a90>

        Calculate and get the inverse transform data for all subswaths in parallel without saving it:

        >>> trans_dat_inv_parallel(interactive=True)
        [<dask.delayed.Delayed at 0x7f8d13a69a90>, <dask.delayed.Delayed at 0x7f8d13a69b70>]
        """
        import dask

        # process all the subswaths
        subswaths = self.get_subswaths()
        delayeds = []
        for subswath in subswaths:
            delayed = self.trans_dat_inv(subswath=subswath, interactive=interactive, **kwargs)
            if not interactive:
                pbar = tqdm_dask(dask.persist(delayed), desc=f'Radar Inverse Transform Computing sw{subswath}')
            else:
                delayeds.append(delayed)

        if interactive:
            return delayeds[0] if len(delayeds)==1 else delayeds

        # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
        import gc; gc.collect()