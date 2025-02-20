# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_export import Stack_export

class Stack(Stack_export):

    def __repr__(self):
        return 'Object %s %d bursts %d dates' % (self.__class__.__name__, len(self.ds), len(self.ds[0].date))

    def to_dataset(self):
        import numpy as np
        import xarray as xr
        data = xr.concat(xr.align(*self.ds, join='outer'), dim='stack_dim').mean('stack_dim')
        return data

    def __init__(self, basedir, pattern_burst='*_*_?W?', pattern_date = '????????.nc', scale = 2.5e-07):
        """
        Initialize an instance of the Stack class.
        """
        import numpy as np
        import xarray as xr
        import glob
        import os

        self.basedir = basedir
        
        bursts = glob.glob(pattern_burst, root_dir=self.basedir)
        datas = []
        for burst in bursts:
            data = xr.open_mfdataset(
                os.path.join(self.basedir, burst, pattern_date),
                engine=self.netcdf_engine_read,
                format=self.netcdf_format,
                parallel=True,
                concat_dim='date',
                chunks={'date': 1, 'y': self.chunksize, 'x': self.chunksize},
                combine='nested',
            )
            #print (data)
            # zero in np.int16 type means NODATA
            data = self.spatial_ref(
                scale*(data.re.astype(np.float32) + 1j*data.im.astype(np.float32)).where(data.re != 0).rename('data'),
                data
            )\
            .to_dataset(name='data')\
            .assign({v: data[v] for v in data.data_vars if v not in ['im', 're']})\
            .assign_attrs({'burst': burst})
            datas.append(data)
            del data
    
        self.ds = datas

    def baseline_table(self):
        import xarray as xr
        return xr.concat([ds.BPR for ds in self.ds], dim='burst').mean('burst').to_dataframe()[['BPR']]

    def baseline_pairs(self, days=None, meters=None, invert=False):
        """
        Generates a sorted list of baseline pairs.
        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the sorted list of baseline pairs with reference and repeat dates,
            timelines, and baselines.
    
        """
        import numpy as np
        import pandas as pd
        
        if days is None:
            # use large number for unlimited time interval in days
            days = 1e6
    
        tbl = self.baseline_table()
        data = []
        for line1 in tbl.itertuples():
            counter = 0
            for line2 in tbl.itertuples():
                #print (line1, line2)
                if not (line1.Index < line2.Index and (line2.Index - line1.Index).days < days + 1):
                    continue
                if meters is not None and not (abs(line1.BPR - line2.BPR)< meters + 1):
                    continue
    
                counter += 1
                if not invert:
                    data.append({'ref':line1.Index, 'rep': line2.Index,
                                 'ref_baseline': np.round(line1.BPR, 2),
                                 'rep_baseline': np.round(line2.BPR, 2)})
                else:
                    data.append({'ref':line2.Index, 'rep': line1.Index,
                                 'ref_baseline': np.round(line2.BPR, 2),
                                 'rep_baseline': np.round(line1.BPR, 2)})
    
        df = pd.DataFrame(data).sort_values(['ref', 'rep'])
        return df.assign(pair=[f'{ref} {rep}' for ref, rep in zip(df['ref'].dt.date, df['rep'].dt.date)],
                         baseline=df.rep_baseline - df.ref_baseline,
                         duration=(df['rep'] - df['ref']).dt.days,
                         rel=np.datetime64('nat'))
