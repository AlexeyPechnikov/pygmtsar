#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_sbas_gmtsar import SBAS_sbas_gmtsar
from .PRM import PRM

class SBAS_sbas(SBAS_sbas_gmtsar):

    def baseline_table(self, n_jobs=-1, debug=False):
        import pandas as pd
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        import os

        # use any subswath (the 1st one here) to produce the table
        subswath = self.get_subswaths()[0]
        datetimes = self.df[self.df.subswath==subswath].datetime

        def get_filename(dt):
            _, stem = self.multistem_stem(subswath, dt)
            filename = os.path.join(self.basedir, f'{stem}.PRM')
            return filename
    
        def ondemand(date, dt):
            if not os.path.exists(get_filename(dt)):
                self.make_s1a_tops(subswath, date, debug=debug)

        # generate PRM, LED if needed
        #for (date, dt) in datetimes.iteritems():
        #    #print (dt, date)
        #    ondemand(dt)
        with self.tqdm_joblib(tqdm(desc='PRM generation', total=len(datetimes))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(ondemand)(date, dt) for (date, dt) in datetimes.iteritems())
    
        # after merging use unmerged subswath PRM files
        # calc_dop_orb() required for SAT_baseline
        master_dt = datetimes[self.master]
        prm_ref = PRM().from_file(get_filename(master_dt)).calc_dop_orb(inplace=True)
        data = []
        for (date, dt) in datetimes.iteritems():
            # after merging use unmerged subswath PRM files
            prm_rep = PRM().from_file(get_filename(dt))
            ST0 = prm_rep.get('SC_clock_start')
            DAY = int(ST0 % 1000)
            YR = int(ST0/1000) - 2014
            YDAY = YR * 365 + DAY
            #print (f'YR={YR}, DAY={DAY}')
            BPL, BPR = prm_ref.SAT_baseline(prm_rep).get('B_parallel', 'B_perpendicular')
            data.append({'date':date, 'ST0':ST0, 'YDAY':YDAY, 'BPL':BPL, 'BPR':BPR})
            #print (date, ST0, YDAY, BPL, BPR)
        return pd.DataFrame(data).set_index('date')

    # returns sorted baseline pairs
    def baseline_pairs(self, days=100, meters=150, invert=False, n_jobs=-1, debug=False):
        import numpy as np
        import pandas as pd
     
        tbl = self.baseline_table(n_jobs=n_jobs, debug=debug)
        data = []
        for line1 in tbl.itertuples():
            for line2 in tbl.itertuples():
                #print (line1, line2)
                if not (line1.YDAY < line2.YDAY and line2.YDAY - line1.YDAY < days):
                    continue
                if not (abs(line1.BPR - line2.BPR)< meters):
                    continue

                if not invert:
                    data.append({'ref_date':line1.Index, 'rep_date': line2.Index,
                                 'ref_timeline': np.round(line1.YDAY/365.25+2014, 2), 'ref_baseline': np.round(line1.BPR, 2),
                                 'rep_timeline': np.round(line2.YDAY/365.25+2014, 2), 'rep_baseline': np.round(line2.BPR, 2)})
                else:
                    data.append({'ref_date':line2.Index, 'rep_date': line1.Index,
                                 'ref_timeline': np.round(line2.YDAY/365.25+2014, 2), 'ref_baseline': np.round(line2.BPR, 2),
                                 'rep_timeline': np.round(line1.YDAY/365.25+2014, 2), 'rep_baseline': np.round(line1.BPR, 2)})

        return pd.DataFrame(data).sort_values(['ref_date', 'rep_date'])

    def sbas_parallel(self, pairs=None, mask=None, data='detrend', corr='corr', chunks=None, n_jobs=-1):
        import xarray as xr
        import numpy as np
        import pandas as pd
        from tqdm.auto import tqdm
        import joblib
        import os

        # compress 3d output following the processing blocks
        netcdf_compression = self.compression.copy()
        netcdf_compression['chunksizes'] = (1, self.chunksize, self.chunksize)

        if pairs is None:
            pairs = self.find_pairs()
        elif isinstance(pairs, pd.DataFrame):
            pairs = pairs.values
        else:
            pairs = np.asarray(pairs)
        # define all the dates as unique reference and repeat dates
        dates = np.unique(pairs.flatten())
    
        # source grids lazy loading
        if isinstance(corr, str):
            corr = self.open_grids(pairs, corr, chunks=chunks)
        if isinstance(data, str):
            data = self.open_grids(pairs, data, chunks=chunks)

        # crop correlation grid like to unwrap grid which may be defined smaller
        corr = corr.reindex_like(data)
        
        # mask can be sparse and limit work area
        if mask is not None:
            data = xr.where(mask>0, data.reindex_like(mask), np.nan)
            corr   = xr.where(mask>0, corr.reindex_like(mask),   np.nan)
    
        # here are one row for every interferogram and one column for every date
        matrix = []
        for pair in pairs:
            mrow = [date>pair[0] and date<=pair[1] for date in dates]
            matrix.append(mrow)
        matrix = np.stack(matrix).astype(int)

        # single-pixel processing function
        def fit(x, w, matrix):
            #return np.zeros(5)
            # ignore pixels where correlation is not defined
            if np.any(np.isnan(w)):
                return np.nan * np.zeros(matrix.shape[1])
            # fill nans as zeroes and set corresponding weight to 0
            nanmask = np.where(np.isnan(x))
            if nanmask[0].size > 0:
                # function arguments are read-only
                x = x.copy()
                w = w.copy()
                x[nanmask] = 0
                w[nanmask] = 0
                # check if x has enough valid values
                if x.size - nanmask[0].size < matrix.shape[1]:
                    return np.nan * np.zeros(matrix.shape[1])
            # least squares solution
            W = (w/np.sqrt(1-w**2))
            model = np.linalg.lstsq(matrix * W[:,np.newaxis], x * W, rcond=None)
            #print ('model', model)
            return model[0]

        # xarray wrapper
        models = xr.apply_ufunc(
            fit,
            data.chunk(dict(pair=-1)),
            corr.chunk(dict(pair=-1)),
            input_core_dims=[['pair'],['pair']],
            exclude_dims=set(('pair',)),
            dask='parallelized',
            vectorize=True,
            output_dtypes=[np.float32],
            output_core_dims=[['date']],
            dask_gufunc_kwargs={'output_sizes': {'date': len(dates)}},
            kwargs={'matrix': matrix}
        )
        # define dates axis
        models['date'] = dates
        # set the stack index to be first
        models = models.transpose('date',...)
        # define the output grid filename
        model_filename = self.get_filenames(None, None, 'disp')
        # cleanup
        if os.path.exists(model_filename):
            os.remove(model_filename)
    
        ts, ys, xs = models.data.blocks.shape
        assert ts == 1, 'Date chunks count should be equal to 1'
        tchunks, ychunks, xchunks = models.chunks
        coordt = models.date
        coordy = np.array_split(models.y, np.cumsum(ychunks))
        coordx = np.array_split(models.x, np.cumsum(xchunks))
    
        def func(iy, ix):
            chunk_filename = self.get_filenames(None, None, f'disp_chunk_{iy}_{ix}')
            if os.path.exists(chunk_filename):
                os.remove(chunk_filename)
            # lazy dask array
            data = models.data.blocks[0,iy,ix]
            # wrap the array
            da = xr.DataArray(data,
                              dims=['date','y','x'],
                              coords={'date': coordt, 'y':coordy[iy], 'x':coordx[ix]})\
                 .rename('displacement')
            # compute and save to NetCDF using chunks of original coordinates
            da.to_netcdf(chunk_filename,
                         unlimited_dims=['y','x'],
                         encoding={'displacement': netcdf_compression},
                         engine=self.engine,
                         compute=True)
            return chunk_filename
    
        # process all the chunks
        with self.tqdm_joblib(tqdm(desc='Computing', total=ys*xs)) as progress_bar:
            filenames = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(func)(iy, ix) \
                                                     for iy in range(ys) for ix in range(xs))
    
        # rebuild the datasets to user-friendly format
        das = [xr.open_dataarray(f, engine=self.engine, chunks=self.chunksize) for f in filenames]
        if xr.__version__ == '0.19.0':
            # for Google Colab
            das = xr.merge(das)
        else:
            # for modern xarray versions
            das = xr.combine_by_coords(das)

        # add subswath prefix
        subswath = self.get_subswath()

        def output(dt):
            filename = os.path.join(self.basedir, f'F{subswath}_disp_{dt}.grd'.replace('-',''))
            if os.path.exists(filename):
                os.remove(filename)
            das.sel(date=dt).to_netcdf(filename,
                        encoding={'displacement': self.compression},
                        engine=self.engine)

        # saving all the grids
        with self.tqdm_joblib(tqdm(desc='Saving', total=len(das.date))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(output)(dt) for dt in das.date.values)

        # cleanup
        for filename in filenames:
            if os.path.exists(filename):
                os.remove(filename)

    #intf.tab format:   unwrap.grd  corr.grd  ref_id  rep_id  B_perp 
    def intftab(self, baseline_pairs):
        import numpy as np
        import datetime

        # return a single subswath for None
        subswath = self.get_subswath()

        outs = []
        for line in baseline_pairs.itertuples():
            #print (line)
            ref = line.ref_date.replace('-','')
            jref = datetime.datetime.strptime(line.ref_date, '%Y-%m-%d').strftime('%Y%j')
            rep = line.rep_date.replace('-','')
            jrep = datetime.datetime.strptime(line.rep_date, '%Y-%m-%d').strftime('%Y%j')
            bperp = np.round(line.rep_baseline - line.ref_baseline, 2)
            outs.append(f'F{subswath}_{ref}_{rep}_unwrap.grd F{subswath}_{ref}_{rep}_corr.grd {jref} {jrep} {bperp}')
        return '\n'.join(outs) + '\n'

    def scenetab(self, baseline_pairs):
        import numpy as np
        import datetime

        mst = datetime.datetime.strptime(self.master, '%Y-%m-%d').strftime('%Y%j')
        #print (self.master, mst)
        outs = []
        for line in baseline_pairs.itertuples():
            #print (line)
            ref = datetime.datetime.strptime(line.ref_date, '%Y-%m-%d').strftime('%Y%j')
            yday_ref = np.round((line.ref_timeline - 2014)*365.25)
            outs.append(f'{ref} {yday_ref}')
            rep = datetime.datetime.strptime(line.rep_date, '%Y-%m-%d').strftime('%Y%j')
            yday_rep = np.round((line.rep_timeline - 2014)*365.25)
            outs.append(f'{rep} {yday_rep}')
        outs = np.unique(outs)
        return '\n'.join([out for out in outs if out.split(' ')[0]==mst]) + '\n' + \
               '\n'.join([out for out in outs if out.split(' ')[0]!=mst]) + '\n'


    @staticmethod
    def triplets2pairs(triplets, pairs):
        import pandas as pd
    
        data = []
        for triplet in triplets.itertuples():
            data.append({'ref_date': triplet.A, 'rep_date': triplet.B})
            data.append({'ref_date': triplet.B, 'rep_date': triplet.C})
            data.append({'ref_date': triplet.A, 'rep_date': triplet.C})
        tripairs = pd.DataFrame(data).sort_values(['ref_date', 'rep_date']).drop_duplicates()
        idx = tripairs.set_index(['ref_date', 'rep_date']).index
        return pairs.set_index(['ref_date', 'rep_date']).loc[idx].reset_index()

    # returns sorted baseline triplets
    @staticmethod
    def pairs2triplets(pairs, invert=False):
        import pandas as pd

        data = []
        pairs_a = pairs
        for line_a in pairs_a.itertuples():
            #print (line_a)
            date_a_ref = line_a.ref_date
            date_a_rep = line_a.rep_date
            pairs_b = pairs[pairs.ref_date==date_a_rep]
            for line_b in pairs_b.itertuples():
                #print (line_b)
                date_b_ref = line_b.ref_date
                date_b_rep = line_b.rep_date
                pairs_c = pairs[(pairs.rep_date==date_b_rep)&(pairs.ref_date==date_a_ref)]
                for line_c in pairs_c.itertuples():
                    #print (line_c)
                    date_c_ref = line_c.ref_date
                    date_c_rep = line_c.rep_date
                    #print (date_a_ref, date_a_rep, date_b_rep)
                    data.append({'A': date_a_ref, 'B': date_a_rep, 'C': date_b_rep})
        return pd.DataFrame(data).sort_values(['A', 'B', 'C'])
