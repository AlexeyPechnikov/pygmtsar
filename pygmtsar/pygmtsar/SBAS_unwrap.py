#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar
from .SBAS_unwrap_snaphu import SBAS_unwrap_snaphu

class SBAS_unwrap(SBAS_unwrap_snaphu):

    def unwrap_parallel(self, pairs, n_jobs=-1, **kwargs):
        import pandas as pd
        from tqdm.auto import tqdm
        import joblib
        import os

        def unwrap(subswath, pair, **kwargs):
            # define unique tiledir name for parallel processing
            if 'conf' in kwargs:
                dirname = f'F{subswath}_{"_".join(pair).replace("-","")}_snaphu_tiledir'
                dirpath = os.path.join(self.basedir, dirname)
                kwargs['conf'] += f'    TILEDIR {dirpath}'
            return self.unwrap(subswath, pair, **kwargs)

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        subswaths = self.get_subswaths()

        with self.tqdm_joblib(tqdm(desc='Unwrapping', total=len(pairs)*len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(unwrap)(subswath, pair, interactive=False, **kwargs) \
                                           for subswath in subswaths for pair in pairs)
