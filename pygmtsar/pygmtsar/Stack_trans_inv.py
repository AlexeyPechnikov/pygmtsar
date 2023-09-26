# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_trans import Stack_trans
from .tqdm_dask import tqdm_dask

class Stack_trans_inv(Stack_trans):

    def get_trans_inv(self):
        """
        Retrieve the inverse transform data.

        This function opens a NetCDF dataset, which contains data mapping from radar
        coordinates to geographical coordinates (from azimuth-range to latitude-longitude domain).

        Parameters
        ----------
        Returns
        -------
        xarray.Dataset
            An xarray dataset with the transform data.

        Examples
        --------
        Get the inverse transform data:
        get_trans_inv()
        """
        return self.open_grid('trans_inv')

    def trans_inv(self, coarsen=1, method='linear', interactive=False):
        """
        Retrieve or calculate the transform data. This transform data is then saved as
            a NetCDF file for future use.

            This function generates data mapping from radar coordinates to geographical coordinates.
            The function uses the direct transform data.

        Parameters
        ----------
        coarsen(jdec, idec) : (int, int) , optional
            The decimation factor in the azimuth and range direction. Default is 2.
        interactive : bool, optional
            If True, the computation will be performed interactively and the result will be returned as a delayed object.
            Default is False.
        """
        import dask
        import xarray as xr
        import numpy as np
        #import os

        # expand simplified definition
        if not isinstance(coarsen, (list,tuple, np.ndarray)):
            coarsen = (coarsen, coarsen)

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

            block_mask = ((trans_dat.amin<=azis_max)&(trans_dat.amax>=azis_min)&(trans_dat.rmin<=rngs_max)&(trans_dat.rmax>=rngs_min)).values
            blocks_ys, blocks_xs = np.meshgrid(trans_dat.block_azi, trans_dat.block_rng, indexing='ij')

            # process chunks
            block_azis = []
            block_rngs = []
            block_lts  = []
            block_lls  = []
            block_eles = []
            for block_y, block_x in zip(blocks_ys[block_mask], blocks_xs[block_mask]):
                # define coordinates
                block_lt, block_ll = [block.ravel() for block in np.meshgrid(lt_blocks[block_y], ll_blocks[block_x], indexing='ij')]
                # extract variables
                block_azi = trans_dat['azi'].data.blocks[block_y, block_x].compute(n_workers=1).ravel()
                block_rng = trans_dat['rng'].data.blocks[block_y, block_x].compute(n_workers=1).ravel()
                block_ele = trans_dat['ele'].data.blocks[block_y, block_x].compute(n_workers=1).ravel()
                # mask valid values
                mask = (block_azi>=azis_min)&(block_azi<=azis_max)&(block_rng>=rngs_min)&(block_rng<=rngs_max)
                block_lts.append(block_lt[mask])
                block_lls.append(block_ll[mask])
                block_azis.append(block_azi[mask])
                block_rngs.append(block_rng[mask])
                block_eles.append(block_ele[mask])
                # cleanup
                del block_lt, block_ll, block_azi, block_rng, block_ele, mask
            # merge extracted results
            return (np.concatenate(block_azis),
                    np.concatenate(block_rngs),
                    np.concatenate(block_lts),
                    np.concatenate(block_lls),
                    np.concatenate(block_eles)
                   )
        
        # griddata interpolation is easy and provides multiple methods
        @dask.delayed
        def trans_dat_inv_block(data, index, azis, rngs):
            from scipy.interpolate import griddata

            #block_azi, block_rng, block_lt, block_ll, block_ele = data
            # coordinates azi, rng
            points = np.column_stack([data[0], data[1]])
            grid = griddata(points, data[index], (azis[None, :], rngs[:, None]), method='linear').astype(np.float32)
            output = grid.reshape((rngs.size, azis.size)).T
            del points, grid
            return output

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

    #     # griddata interpolation is easy and provides multiple methods
    #     @dask.delayed
    #     def trans_dat_inv_block(data, azis, rngs):
    #         from scipy.interpolate import griddata

    #         block_azi, block_rng, block_lt, block_ll, block_ele = data
    #         points = np.column_stack([block_azi, block_rng])

    #         output = []
    #         for block in [block_lt, block_ll, block_ele]:
    #             grid = griddata(points, block, (azis[None, :], rngs[:, None]), method='linear').astype(np.float32)
    #             output.append(grid.reshape((rngs.size, azis.size)).T)
    #             del grid

    #         del block_ele, block_azi, block_rng, block_lt, block_ll, points
    #         return np.asarray(output)

        # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        trans_dat = self.get_trans()
        # check that the dataset is opened with the same chunksize as it is used for creation
        assert trans_dat.azi.data.numblocks == trans_dat.amax.shape, 'Chunksize is differ to the original for trans_dat'
        # materialize indices
        trans_dat['amin'] = trans_dat.amin.compute()
        trans_dat['amax'] = trans_dat.amax.compute()
        trans_dat['rmin'] = trans_dat.rmin.compute()
        trans_dat['rmax'] = trans_dat.rmax.compute()
        # convert axes to Dask arrays with the same chunks
        chunks = trans_dat.azi.data.chunks
        lt_blocks = np.array_split(trans_dat['lat'].values, np.cumsum(chunks[0])[:-1])
        ll_blocks = np.array_split(trans_dat['lon'].values, np.cumsum(chunks[1])[:-1])

        # define topo_ra grid
        azis, rngs = self.define_trans_grid(coarsen)

        # build topo_ra grid by chunks
        # split to equal chunks and rest
        azis_blocks = np.array_split(azis, np.arange(0, azis.size, self.chunksize)[1:])
        rngs_blocks = np.array_split(rngs, np.arange(0, rngs.size, self.chunksize)[1:])
        #print ('azis_blocks.size', len(azis_blocks), 'rngs_blocks.size', len(rngs_blocks))
        #print ('lats_blocks[0]', lats_blocks[0])

        block_lts_total  = []
        block_lls_total  = []
        block_eles_total = []
        for rngs_block in rngs_blocks:
            block_lts  = []
            block_lls  = []
            block_eles = []
            for azis_block in azis_blocks:
                data = trans_dat_inv_block_prepare(azis_block, rngs_block, self.chunksize)
                block_lt  = dask.array.from_delayed(trans_dat_inv_block(data, 2, azis_block, rngs_block),
                                                shape=(azis_block.size, rngs_block.size), dtype=np.float32)
                block_lts.append(block_lt.transpose(1,0))
                block_ll  = dask.array.from_delayed(trans_dat_inv_block(data, 3, azis_block, rngs_block),
                                                shape=(azis_block.size, rngs_block.size), dtype=np.float32)
                block_lls.append(block_ll.transpose(1,0))
                block_ele = dask.array.from_delayed(trans_dat_inv_block(data, 4, azis_block, rngs_block),
                                                shape=(azis_block.size, rngs_block.size), dtype=np.float32)
                block_eles.append(block_ele.transpose(1,0))
                del block_lt, block_ll, block_ele
            block_lts_total.append(block_lts)
            block_lls_total.append(block_lls)
            block_eles_total.append(block_eles)
            del block_lts, block_lls, block_eles

        coords = {'y': azis, 'x': rngs}
        lt  = xr.DataArray(dask.array.block(block_lts_total).transpose(1,0),  coords=coords)
        ll  = xr.DataArray(dask.array.block(block_lls_total).transpose(1,0),  coords=coords)
        ele = xr.DataArray(dask.array.block(block_eles_total).transpose(1,0), coords=coords)
        del block_lts_total, block_lls_total, block_eles_total

        trans_inv = xr.Dataset({'lt': lt, 'll': ll, 'ele': ele})
        del lt, ll, ele

        # calculate indices
        trans_inv['lt_min'] = trans_inv.lt.min('x')
        trans_inv['lt_max'] = trans_inv.lt.max('x')
        trans_inv['ll_min'] = trans_inv.ll.min('y')
        trans_inv['ll_max'] = trans_inv.ll.max('y')

        if interactive:
            return trans_inv
        # rename to save lazy NetCDF preventing broken coordinates (y,y)
        return self.save_grid(trans_inv.rename({'y': 'a', 'x': 'r'}),
                              'trans_inv', f'Radar Inverse Transform Computing')
