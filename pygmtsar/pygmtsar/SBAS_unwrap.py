#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar
from .SBAS_unwrap_snaphu import SBAS_unwrap_snaphu

class SBAS_unwrap(SBAS_unwrap_snaphu):

    def unwrap_parallel(self, pairs=None, mask=None, n_jobs=-1, **kwargs):
        import xarray as xr
        import pandas as pd
        from tqdm.auto import tqdm
        import joblib
        import os

        # for now (Python 3.10.10 on MacOS) joblib loads the code from disk instead of copying it
        kwargs['chunksize'] = self.chunksize

        def unwrap_tiledir(pair, **kwargs):
            # define unique tiledir name for parallel processing
            if 'conf' in kwargs:
                dirpath = self.get_filenames(None, [pair], 'snaphu_tiledir')[0][:-4]
                kwargs['conf'] += f'    TILEDIR {dirpath}'
            return self.unwrap(pair, **kwargs)

        if pairs is None:
            pairs = self.find_pairs()
        elif isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        # materialize lazy mask
        if mask is not None and isinstance(mask, xr.DataArray):
            mask_filename = self.get_filenames(None, None, 'unwrapmask')
            if os.path.exists(mask_filename):
                os.remove(mask_filename)
            # workaround to save NetCDF file correct
            mask.rename('mask').rename({'y':'a','x':'r'}).\
                to_netcdf(mask_filename, encoding={'mask': self.compression(chunksize=self.chunksize)}, engine=self.engine)
            kwargs['mask'] = 'unwrapmask'

        # save results to NetCDF files
        kwargs['interactive'] = False

        with self.tqdm_joblib(tqdm(desc='Unwrapping', total=len(pairs))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(unwrap_tiledir)(pair, **kwargs) for pair in pairs)
