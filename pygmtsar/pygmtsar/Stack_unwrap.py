# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_unwrap_snaphu import Stack_unwrap_snaphu
# required for function decorators
from numba import jit
# import directive is not compatible to numba
import numpy as np

class Stack_unwrap(Stack_unwrap_snaphu):

    @staticmethod
    def wrap(data_pairs):
        import xarray as xr
        import numpy as np
        import dask

        if isinstance(data_pairs, xr.DataArray):
            return xr.DataArray(dask.array.mod(data_pairs.data + np.pi, 2 * np.pi) - np.pi, data_pairs.coords)\
                .rename(data_pairs.name)
        return np.mod(data_pairs + np.pi, 2 * np.pi) - np.pi

    # weight is used only to filter NaNs in data
    @staticmethod
    @jit(nopython=True, nogil=True)
    def unwrap_pairs(data, weight, matrix):
        # import directive is not compatible to numba
        #import numpy as np

        # initial data size
        data_size = data.size

        if np.all(np.isnan(data)) or (weight is not None and np.all(np.isnan(weight))):
            return np.nan * np.zeros(data_size)

        # the data variable will be modified and returned as the function output
        data = data.copy()
        if weight is not None:
            nanmask = np.isnan(data) | np.isnan(weight)
        else:
            nanmask = np.isnan(data)
        if np.all(nanmask):
            return np.nan * np.zeros(data_size)
        data = data[~nanmask].copy()
        if weight is not None:
            weight = weight[~nanmask]
        matrix = matrix[~nanmask,:]

        #data = phase_pairs.values
        pair_sum = matrix.sum(axis=1)
        # tidal displacements for SBAS pairs
        #tidal_pairs = (matrix * np.concatenate([[0], np.diff(tidal)])).sum(axis=1)
        # wrapped phase pairs (phase_tidal)
        #data = (np.mod((phase_pairs - tidal_pairs) + np.pi, 2 * np.pi) - np.pi)

        # pairs do not require unwrapping
        pairs_ok = []

        # check all compound pairs vs single pairs: only detect all not wrapped
        for ndate in np.unique(pair_sum)[1:]:
            for pair_idx in np.where(pair_sum==ndate)[0]:
                #print (pair_idx, matrix[pair_idx])
                matching_columns = matrix[pair_idx] == 1
                #print ('matching_columns', matching_columns)
                # works for dates = 2
                #matching_rows = ((matrix[:, matching_columns] == 1).sum(axis=1) == 1)&((matrix[:, ~matching_columns] == 1).sum(axis=1) == 0)
                matching_rows = ((matrix[:, matching_columns] >= 1).sum(axis=1) == 1)&((matrix[:, ~matching_columns] == 1).sum(axis=1) == 0)
                #print ('matching_rows', matching_rows)
                matching_matrix = matrix[:,matching_columns] * matching_rows[:,None]
                #row_indices = np.where(matching_rows)[0]
                #print ('row_indices', row_indices)
                #with np.printoptions(threshold=np.inf, linewidth=np.inf):
                #    print (matching_matrix.T)
                value = data[pair_idx]
                values = (matching_matrix * data[:,None])
                #print ('value', value, '?=', values.sum())    
                jump = int(np.round((values.sum() - value) / (2*np.pi)))
                if jump == 0:
                    pairs_ok.append(pair_idx)
                    pairs_ok.extend(np.where(matching_rows)[0])
                    #print ()
                    continue

                #print ('jump', jump)
                if abs(value) > abs(value + 2*np.pi*jump):
                    #print (f'JUMP {ndate}:', jump)
                    #jump, phase.values[pair_idx] + 2*np.pi*jump
                    #print ('value', value, '=>', value + 2*np.pi*jump, '=', values.sum())
                    pass
                else:
                    #print (f'JUMP singles:', jump)
                    pass
                #print (values[matching_rows])
                # check for wrapping
                valid_values = values[matching_rows].ravel()
                #maxdiff = abs(np.diff(valid_values[valid_values!=0])).max()
                #print ('maxdiff', maxdiff, maxdiff >= np.pi)
                #print ()
                #break
        #print ('pairs_ok', pairs_ok)
        #print ('==================')

        # check all compound pairs vs single pairs: fix wrapped compound using not wrapped singles only
        for ndate in np.unique(pair_sum)[1:]:
            for pair_idx in np.where(pair_sum==ndate)[0]:
                if pair_idx in pairs_ok:
                    continue
                #print (pair_idx, matrix[pair_idx])
                matching_columns = matrix[pair_idx] == 1
                #print ('matching_columns', matching_columns)
                # works for dates = 2
                #matching_rows = ((matrix[:, matching_columns] == 1).sum(axis=1) == 1)&((matrix[:, ~matching_columns] == 1).sum(axis=1) == 0)
                matching_rows = ((matrix[:, matching_columns] >= 1).sum(axis=1) == 1)&((matrix[:, ~matching_columns] == 1).sum(axis=1) == 0)
                #print ('matching_rows', matching_rows)
                matching_matrix = matrix[:,matching_columns] * matching_rows[:,None]
                #row_indices = np.where(matching_rows)[0]
                #print ('row_indices', row_indices)
                #with np.printoptions(threshold=np.inf, linewidth=np.inf):
                #    print (matching_matrix.T)
                #print ('matching_matrix', matching_matrix.shape)

                pairs_single_not_valid = [pair for pair in np.where(matching_rows)[0] if not pair in pairs_ok]
                if len(pairs_single_not_valid) > 0:
                    # some of single-pairs requires unwrapping, miss the compound segment processing
                    #print ('ERROR')
                    #print ()
                    continue

                value = data[pair_idx]
                values = (matching_matrix * data[:,None])
                #print ('value', value, '?=', values.sum())    
                jump = int(np.round((values.sum() - value) / (2*np.pi)))
                #print (f'JUMP {ndate}:', jump)
                data[pair_idx] += 2*np.pi*jump
                pairs_ok.append(pair_idx)
                #print ()
                continue
        #print ('==================')

        # check all compound pairs vs single pairs
        for ndate in np.unique(pair_sum)[1:]:
            for pair_idx in np.where(pair_sum==ndate)[0]:
                if pair_idx in pairs_ok:
                    continue
                #print (pair_idx, matrix[pair_idx])
                matching_columns = matrix[pair_idx] == 1
                #print ('matching_columns', matching_columns)
                # works for dates = 2
                #matching_rows = ((matrix[:, matching_columns] == 1).sum(axis=1) == 1)&((matrix[:, ~matching_columns] == 1).sum(axis=1) == 0)
                matching_rows = ((matrix[:, matching_columns] >= 1).sum(axis=1) == 1)&((matrix[:, ~matching_columns] == 1).sum(axis=1) == 0)
                #print ('matching_rows', matching_rows)
                matching_matrix = matrix[:,matching_columns] * matching_rows[:,None]
                #row_indices = np.where(matching_rows)[0]
                #print ('row_indices', row_indices)
                #with np.printoptions(threshold=np.inf, linewidth=np.inf):
                #    print (matching_matrix.T)
                value = data[pair_idx]
                values = (matching_matrix * data[:,None])
                #print ('value', value, '???=', values.sum())    
                jump = int(np.round((values.sum() - value) / (2*np.pi)))
                #print ('jump', jump)
                if jump != 0:
                    #print (f'JUMP {ndate}:', jump)
                    #jump, phase.values[pair_idx] + 2*np.pi*jump
                    #print ('value', value, '=>', value + 2*np.pi*jump, '=', values.sum())
                    pass
                #print (values[matching_rows])
                # unwrap
                data[pair_idx] += 2*np.pi*jump
                pairs_ok.append(pair_idx)
                pairs_ok.extend(np.where(matching_rows)[0])
                #print ()
                #break

        # validity mask
        mask = [idx in pairs_ok for idx in range(data.size)]
        data = np.where(mask, data, np.nan)
        out = np.nan * np.zeros(data_size)
        out[~nanmask] = data
        return out

    def unwrap1d(self, data, weight=None):
        import xarray as xr
        import numpy as np

        chunks_z, chunks_y, chunks_x = data.chunks
        if np.max(chunks_y) > self.netcdf_chunksize or np.max(chunks_x) > self.netcdf_chunksize:
            print (f'Note: data chunk size ({np.max(chunks_y)}, {np.max(chunks_x)}) is too large for stack processing')
            chunks_y = chunks_x = self.netcdf_chunksize//2
            print (f'Note: auto tune data chunk size to a half of NetCDF chunk: ({chunks_y}, {chunks_x})')
        else:
            # use the existing data chunks size
            chunks_y = np.max(chunks_y)
            chunks_x = np.max(chunks_x)

        pairs = self.get_pairs(data)
        matrix = self.lstsq_matrix(pairs)

        # xarray wrapper
        model = xr.apply_ufunc(
            self.unwrap_pairs,
            data.chunk(dict(pair=-1, y=chunks_y, x=chunks_x)),
            weight.chunk(dict(pair=-1, y=chunks_y, x=chunks_x)) if weight is not None else weight,
            dask='parallelized',
            vectorize=True,
            input_core_dims=[['pair'], ['pair']] if weight is not None else [['pair'], []],
            output_core_dims=[['pair']],
            output_dtypes=[np.float32],
            kwargs={'matrix': matrix}
        ).transpose('pair',...)

        return model.rename('unwrap')

    def unwrap_snaphu(self, phase, weight=None, conf=None, conncomp=False):
        """
        Limit number of processes for tiled multicore SNAPHU configuration:
        with dask.config.set(scheduler='single-threaded'):
            stack.unwrap2d_snaphu()...).phase.compute()
        or
        with dask.config.set(scheduler='threads', num_workers=2):
            stack.unwrap2d_snaphu()...).phase.plot.imshow()
        See for reference: https://docs.dask.org/en/stable/scheduler-overview.html
        """
        import xarray as xr
        import numpy as np
        import dask

        # output dataset variables
        keys = {0: 'phase', 1: 'conncomp'}

        if len(phase.dims) == 2:
            stackvar = None
        else:
            stackvar = phase.dims[0]

        if weight is not None:
            assert phase.shape == weight.shape, 'ERROR: phase and weight variables have different shape'

        def _snaphu(ind):
            ds = self.snaphu(phase.isel({stackvar: ind}) if stackvar is not None else phase,
                             weight.isel({stackvar: ind})  if stackvar is not None and weight is not None else weight,
                             conf=conf, conncomp=conncomp)
            if conncomp:
                # # select the largest connected component
#                 hist = np.unique(ds.conncomp.values, return_counts=True)
#                 idxmax = np.argmax(hist[1])
#                 valmax = hist[0][idxmax]
#                 # select unwrap phase for the largest connected area
#                 return np.stack([ds.phase.where(ds.conncomp==valmax).values, ds.conncomp.values])
                return np.stack([ds.phase.values, ds.conncomp.values])
            return ds.phase.values[None,]

        stack =[]
        for ind in range(len(phase) if stackvar is not None else 1):
            block = dask.array.from_delayed(dask.delayed(_snaphu)(ind),
                        shape=(2 if conncomp else 1, *(phase.shape[1:] if stackvar is not None else phase.shape)),
                        dtype=np.float32)
            stack.append(block)
            del block
        dask_block = dask.array.concatenate(stack)
        del stack
        if stackvar is not None:
            ds = xr.merge([xr.DataArray(dask_block[idx::2 if conncomp else 1], coords=phase.coords).rename(keys[idx])
                               for idx in range(2 if conncomp else 1)])
        else:
            ds = xr.merge([xr.DataArray(block, coords=phase.coords).rename(keys[idx])
                           for idx, block in enumerate(dask_block)])
        del dask_block
        return ds

    def interpolate_nearest(self, data, search_radius_pixels=None):
        """
        Perform nearest neighbor interpolation on each 2D grid in a 3D grid stack.

        Parameters
        ----------
        data : xarray.DataArray
            The input 3D grid stack to be interpolated, with dimensions (pair, y, x).
        search_radius_pixels : int, optional
            The interpolation distance in pixels. If not provided, the default is set to the chunksize of the Stack object.

        Returns
        -------
        xarray.DataArray
            The interpolated 3D grid stack.

        Examples
        --------
        Fill gaps in the specified grid stack using nearest neighbor interpolation:
        stack.interpolate_nearest(intf)

        Notes
        -----
        This method performs nearest neighbor interpolation on each 2D grid (y, x) in a 3D grid stack (pair, y, x). It replaces the NaN values in each 2D grid with the nearest non-NaN values. The interpolation is performed within a specified search radius in pixels for each grid. If a search radius is not provided, the default search radius is set to the chunksize of the Stack object.
        """
        import xarray as xr

        assert data.dims == ('pair', 'y', 'x'), 'Input data must have dimensions (pair, y, x)'

        interpolated = [self.nearest_grid(data.sel(pair=pair), search_radius_pixels) for pair in data.pair]
        return xr.concat(interpolated, dim='pair')

    @staticmethod
    def conncomp_main(data):
        import xarray as xr
        import numpy as np

        # Function to find the mode (most frequent value)
        def find_mode(array):
            values, counts = np.unique(array[~np.isnan(array)], return_counts=True)
            max_count_index = np.argmax(counts)
            return values[max_count_index] if counts.size > 0 else np.nan

        # Apply the function along the 'pair' dimension
        maincomps =  xr.apply_ufunc(find_mode, data.conncomp, input_core_dims=[['y', 'x']],
                                    vectorize=True, dask='parallelized', output_dtypes=[float])
        return data.where(data.conncomp==maincomps)

    def plot_conncomps(self, data, caption='Connected Components', cols=4, size=4, nbins=5, aspect=1.2, y=1.05,
                       vmin=0, vmax=10, cmap='tab10_r'):
        import matplotlib.pyplot as plt

        # multi-plots ineffective for linked lazy data
        fg = data.plot.imshow(
            col='pair',
            col_wrap=cols, size=size, aspect=aspect,
            cmap=cmap, vmin=0, vmax=10
        )
        fg.set_axis_labels('Range', 'Azimuth')
        fg.set_ticks(max_xticks=nbins, max_yticks=nbins)
        fg.fig.suptitle(caption, y=y)
