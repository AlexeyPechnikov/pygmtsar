#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar
from .SBAS_unwrap_snaphu import SBAS_unwrap_snaphu

class SBAS_unwrap(SBAS_unwrap_snaphu):

    def unwrap_parallel(self, pairs=None, n_jobs=-1, **kwargs):
        import pandas as pd
        from tqdm.auto import tqdm
        import joblib
        import os

        def unwrap(pair, **kwargs):
            # define unique tiledir name for parallel processing
            if 'conf' in kwargs:
                dirpath = self.get_filenames(None, [pair], 'snaphu_tiledir')[0][:-4]
                kwargs['conf'] += f'    TILEDIR {dirpath}'
            return self.unwrap(pair, **kwargs)

        if pairs is None:
            pairs = self.find_pairs()
        elif isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        with self.tqdm_joblib(tqdm(desc='Unwrapping', total=len(pairs))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(unwrap)(pair, interactive=False, **kwargs) \
                                           for pair in pairs)
