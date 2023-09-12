# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2022, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_unwrap_snaphu import Stack_unwrap_snaphu

class Stack_unwrap(Stack_unwrap_snaphu):

    @staticmethod
    def wrap(data_pairs):
        import numpy as np

        return np.mod(data_pairs + np.pi, 2 * np.pi) - np.pi

    @staticmethod
    def unwrap1d(data_pairs, matrix):
        import numpy as np

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
                jump = int(((values.sum() - value) / (2*np.pi)).round())
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
                jump = int(((values.sum() - value) / (2*np.pi)).round())
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
                jump = int(((values.sum() - value) / (2*np.pi)).round())
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

    def stack_unwrap1d(self, data):
        import xarray as xr
        import numpy as np
        
        baseline_pairs = self.get_pairs(data)
        matrix = self.lstsq_matrix(baseline_pairs)

        # xarray wrapper
        model = xr.apply_ufunc(
            self.unwrap1d,
            data.chunk(dict(pair=-1)),
            dask='parallelized',
            vectorize=True,
            input_core_dims=[['pair']],
            output_core_dims=[['pair']],
            output_dtypes=[np.float32],
            kwargs={'matrix': matrix}
        ).transpose('pair',...)

        return model.rename('unwrap')

    # TODO
    def stack_snaphu(self):
        return
# 
#     def snaphu_parallel(self, pairs=None, mask=None, chunksize=None, n_jobs=-1, **kwargs):
#         """
#         Unwraps phase using SNAPHU in parallel for multiple interferogram pairs.
# 
#         This function unwraps the phase of multiple interferograms using the Statistical-cost, Network-flow Algorithm
#         for Phase Unwrapping (SNAPHU) with user-defined parameters. The unwrapped phase is saved as grid files
#         in the working directory. This function is designed to process multiple interferograms in parallel, 
#         leveraging the available computational resources.
# 
#         Parameters
#         ----------
#         pairs : list, tuple, array or pandas.DataFrame, optional
#             A list or array of pairs of reference and repeat dates, or a DataFrame with 'ref' and 'rep' columns.
# 
#         mask : str or xarray.DataArray, optional
#             A user-defined mask for radar coordinates, default is None.
# 
#         n_jobs : int, optional
#             The number of jobs to run in parallel. -1 means using all available processors, default is -1.
# 
#         **kwargs : dict
#             Additional keyword arguments to be passed to the unwrap function.
# 
#         Returns
#         -------
#         None
#             This function saves the unwrapped phase to disk as grid files and does not return any values.
# 
#         Examples
#         --------
#         Simplest unwrapping:
#         stack.unwrap_parallel(pairs)
# 
#         Filter low-coherence areas with common correlation threshold 0.75:
#         stack.unwrap_parallel(pairs, threshold=0.075)
# 
#         Unwrap with coherence threshold 0.075 and fill NODATA gaps:
#         interpolator = lambda corr, unwrap: stack.nearest_grid(unwrap).where(corr>=0)
#         stack.unwrap_parallel(pairs, threshold=0.075, func=interpolator)
# 
#         Unwrap with coherence threshold 0.075 and apply land mask:
#         cleaner = lambda corr, unwrap: xr.where(corr>=0.075, unwrap, np.nan)
#         stack.unwrap_parallel(pairs, threshold=0.075, mask=landmask_ra, func=cleaner)
# 
#         Unwrap with coherence threshold 0.075 and use SNAPHU tiling for faster processing and smaller RAM usage:
#         cleaner = lambda corr, unwrap: xr.where(corr>=0.075, unwrap, np.nan)
#         conf = stack.PRM().snaphu_config(NTILEROW=1, NTILECOL=2, ROWOVRLP=200, COLOVRLP=200)
#         stack.unwrap_parallel(pairs, n_jobs=1, threshold=0.075, func=cleaner, conf=conf)
# 
#         Notes
#         -----
#         This function unwraps phase using SNAPHU in parallel for multiple interferogram pairs.
#         The unwrapped phase is saved as grid files in the working directory.
#         Additional keyword arguments can be passed to customize the unwrapping process.
#         """
#         import xarray as xr
#         import pandas as pd
#         from tqdm.auto import tqdm
#         import joblib
#         import os
# 
#         # for now (Python 3.10.10 on MacOS) joblib loads the code from disk instead of copying it
#         kwargs['chunksize'] = chunksize
# 
#         # convert pairs (list, array, dataframe) to 2D numpy array
#         pairs = self.get_pairs(pairs)[['ref', 'rep']].astype(str).values
#         
#         # materialize lazy mask
# #         if mask is not None and isinstance(mask, xr.DataArray):
# #             mask_filename = self.get_filenames(None, 'unwrapmask')
# #             if os.path.exists(mask_filename):
# #                 os.remove(mask_filename)
# #             # workaround to save NetCDF file correct
# #             #if 'y' in mask.dims and 'x' in mask.dims:
# #             #    mask = mask.rename({'y':'a', 'x':'r'})
# #             mask.rename('mask').to_netcdf(mask_filename, encoding={'mask': self.compression(mask.shape, chunksize=chunksize)}, engine=self.engine)
# #             kwargs['mask'] = 'unwrapmask'
# 
#         # save results to NetCDF files
#         kwargs['interactive'] = False
# 
#         with self.tqdm_joblib(tqdm(desc='SNAPHU Unwrapping', total=len(pairs))) as progress_bar:
#             joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.snaphu)(pair, tiledir=f'snaphu_tiledir_{pair[0]}_{pair[1]}', **kwargs) for pair in pairs)
# # 