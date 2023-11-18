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
from numba import jit
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

    @staticmethod
    @jit(nopython=True)
    def unwrap_pairs(data_pairs, matrix):
        # import directive is not compatible to numba
        #import numpy as np

        # NaN values are not permitted; TBD: its possible to exclude some NaNs
        if np.any(np.isnan(data_pairs)):
            return np.nan * np.zeros(data_pairs.shape)

        # the data variable will be modified and returned as the function output
        data = data_pairs.copy()
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
        mask = [idx in pairs_ok for idx in range(data_pairs.size)]
        return np.where(mask, data, np.nan)

    def unwrap1d(self, data):
        import xarray as xr
        import numpy as np
        
        baseline_pairs = self.get_pairs(data)
        matrix = self.lstsq_matrix(baseline_pairs)

        # xarray wrapper
        model = xr.apply_ufunc(
            self.unwrap_pairs,
            data.chunk(dict(pair=-1)),
            dask='parallelized',
            vectorize=True,
            input_core_dims=[['pair']],
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
