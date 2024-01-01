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
        return self.open_cube('trans_inv')

    def compute_trans_inv(self, coarsen, trans='auto', interactive=False):
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

        Note
        ----
        This function operates on the 'trans' grid using NetCDF chunks (specified by 'netcdf_chunksize') rather than
        larger processing chunks. This approach is effective due to on-the-fly index creation for the NetCDF chunks.

        """
        import dask
        import xarray as xr
        import numpy as np
        import warnings
        warnings.filterwarnings('ignore')

        # convert meters or pixels to radar pixels
        coarsen = self.get_coarsen(coarsen)
        # define maximum search radius, radar pixels
        tolerance = 2 * max(coarsen)

        def trans_inv_block(azis, rngs, tolerance, chunksize):
            from scipy.spatial import cKDTree
            # disable "distributed.utils_perf - WARNING - full garbage collections ..."
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()
            import warnings
            warnings.filterwarnings('ignore')

            # required one delta around for nearest interpolation and two for linear
            dazi = np.diff(azis)[0]
            drng = np.diff(rngs)[0]
            azis_min = azis.min() - dazi
            azis_max = azis.max() + dazi
            rngs_min = rngs.min() - drng
            rngs_max = rngs.max() + drng
            del dazi, drng
            #print ('azis_min', azis_min, 'azis_max', azis_max)

            # define valid coordinate blocks 
            block_mask = ((trans_amin<=azis_max)&(trans_amax>=azis_min)&(trans_rmin<=rngs_max)&(trans_rmax>=rngs_min)).values
            block_azi, block_rng = trans_amin.shape
            blocks_ys, blocks_xs = np.meshgrid(range(block_azi), range(block_rng), indexing='ij')
            #assert 0, f'blocks_ys, blocks_xs: {blocks_ys[block_mask]}, {blocks_xs[block_mask]}'
            # extract valid coordinates from the defined blocks
            blocks_trans = []
            blocks_lt = []
            blocks_ll = []
            for block_y, block_x in zip(blocks_ys[block_mask], blocks_xs[block_mask]):
                # coordinates
                block_lt, block_ll = [block.ravel() for block in np.meshgrid(lt_blocks[block_y], ll_blocks[block_x], indexing='ij')]
                # variables
                block_trans = trans.isel(lat=slice(chunksize*block_y,chunksize*(block_y+1)),
                                         lon=slice(chunksize*block_x,chunksize*(block_x+1)))[['azi', 'rng', 'ele']]\
                                   .compute(n_workers=1).to_array().values.reshape(3,-1)
                # select valuable coordinates only
                mask = (block_trans[0,:]>=azis_min)&(block_trans[0,:]<=azis_max)&\
                       (block_trans[1,:]>=rngs_min)&(block_trans[1,:]<=rngs_max)
                # ignore block without valid pixels
                if mask[mask].size > 0:
                    # append valid pixels to accumulators
                    blocks_lt.append(block_lt[mask])
                    blocks_ll.append(block_ll[mask])
                    blocks_trans.append(block_trans[:,mask])
                del block_lt, block_ll, block_trans, mask
            del block_mask, block_azi, block_rng, blocks_ys, blocks_xs

            if len(blocks_lt) == 0:
                # this case is possible when DEM is incomplete, and it is not an error
                return np.nan * np.zeros((3, azis.size, rngs.size), np.float32)

            # TEST
            #return np.nan * np.zeros((3, azis.size, rngs.size), np.float32)

            # valid coordinates
            block_lt = np.concatenate(blocks_lt)
            block_ll = np.concatenate(blocks_ll)
            block_trans = np.concatenate(blocks_trans, axis=1)
            del blocks_lt, blocks_ll, blocks_trans

            # perform index search on radar coordinate grid for the nearest geographic coordinates grid pixel
            grid_azi, grid_rng = np.meshgrid(azis, rngs, indexing='ij')
            tree = cKDTree(np.column_stack([block_trans[0], block_trans[1]]), compact_nodes=False, balanced_tree=False)
            distances, indices = tree.query(np.column_stack([grid_azi.ravel(), grid_rng.ravel()]), k=1, workers=1)
            del grid_azi, grid_rng, tree, cKDTree

            # take the nearest pixels coordinates and elevation
            # the only one index search is required to define all the output variables
            grid_lt = block_lt[indices]
            grid_lt[distances>tolerance] = np.nan
            del block_lt
            grid_ll = block_ll[indices]
            grid_ll[distances>tolerance] = np.nan
            del block_ll
            grid_ele = block_trans[2][indices]
            grid_ele[distances>tolerance] = np.nan
            #print ('distance range', distances.min().round(2), distances.max().round(2))
            #assert distances.max() < 2, f'Unexpectedly large distance between radar and geographic coordinate grid pixels (>=2): {distances.max()}'
            del block_trans, indices, distances

            # pack all the outputs into one 3D array
            return np.asarray([grid_lt, grid_ll, grid_ele]).reshape((3, azis.size, rngs.size))

        if isinstance(trans, str) and trans == 'auto':
            # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
            trans = self.get_trans()

        # calculate indices on the fly
        trans_blocks = trans[['azi', 'rng']].coarsen(lat=self.netcdf_chunksize, lon=self.netcdf_chunksize, boundary='pad')
        #block_min, block_max = dask.compute(trans_blocks.min(), trans_blocks.max())
        # materialize with progress bar indication
        tqdm_dask(trans_blocks_persist := dask.persist(trans_blocks.min(), trans_blocks.max()), desc='Radar Transform Indexing')
        # only convert structure
        block_min, block_max = dask.compute(trans_blocks_persist)[0]
        trans_amin = block_min.azi
        trans_amax = block_max.azi
        trans_rmin = block_min.rng
        trans_rmax = block_max.rng
        del trans_blocks, block_min, block_max
        #print ('trans_amin', trans_amin)

        # split geographic coordinate grid to equal chunks and rest
        #chunks = trans.azi.data.chunks
        #lt_blocks = np.array_split(trans['lat'].values, np.cumsum(chunks[0])[:-1])
        #ll_blocks = np.array_split(trans['lon'].values, np.cumsum(chunks[1])[:-1])
        lt_blocks = np.array_split(trans['lat'].values, np.arange(0, trans['lat'].size, self.netcdf_chunksize)[1:])
        ll_blocks = np.array_split(trans['lon'].values, np.arange(0, trans['lon'].size, self.netcdf_chunksize)[1:])

        # split radar coordinate grid to equal chunks and rest
        azis, rngs = self.define_trans_grid(coarsen)
        azis_blocks = np.array_split(azis, np.arange(0, azis.size, self.netcdf_chunksize)[1:])
        rngs_blocks = np.array_split(rngs, np.arange(0, rngs.size, self.netcdf_chunksize)[1:])
        #print ('azis_blocks.size', len(azis_blocks), 'rngs_blocks.size', len(rngs_blocks))

        blocks_total = []
        for azis_block in azis_blocks:
            blocks = []
            for rngs_block in rngs_blocks:
                block = dask.array.from_delayed(dask.delayed(trans_inv_block, traverse=False)
                                               (azis_block, rngs_block, tolerance, self.netcdf_chunksize),
                                               shape=(3, azis_block.size, rngs_block.size), dtype=np.float32)
                blocks.append(block)
                del block
            blocks_total.append(blocks)
            del blocks

        trans_inv_dask = dask.array.block(blocks_total)
        del blocks_total
        coords = {'y': azis, 'x': rngs}
        trans_inv = xr.Dataset({key: xr.DataArray(trans_inv_dask[idx],  coords=coords)
                                for idx, key in enumerate(['lt', 'll', 'ele'])})
        del trans_inv_dask

        # calculate indices
        trans_inv['lt_min'] = trans_inv.lt.min('x')
        trans_inv['lt_max'] = trans_inv.lt.max('x')
        trans_inv['ll_min'] = trans_inv.ll.min('y')
        trans_inv['ll_max'] = trans_inv.ll.max('y')

        # add target geographic coordinate grid for the user defined spacing (coarsen)
        trans_inv['lat'] = trans['lat'].values
        trans_inv['lon'] = trans['lon'].values

        if interactive:
            return trans_inv
        # rename to save lazy NetCDF preventing broken coordinates (y,y)
        return self.save_cube(trans_inv, 'trans_inv', f'Radar Inverse Transform Computing')
