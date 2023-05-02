#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_sbas_gmtsar import SBAS_sbas_gmtsar
from .PRM import PRM
from .tqdm_dask import tqdm_dask

class SBAS_sbas(SBAS_sbas_gmtsar):

    # single-pixel processing function
    # compute least squares when w = None and weighted least squares otherwise
    @staticmethod
    def lstsq(x, w, matrix):
        import numpy as np

        if w is None:
            # allow not weighted Least Squares
            w = 0.5 * np.ones_like(x)
        elif np.any(np.isnan(w)):
            # ignore pixels where correlation is not defined
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
        try:
            if np.all(w == 1):
                # least squares solution
                model = np.linalg.lstsq(matrix, x, rcond=None)
            else:
                # weighted least squares solution
                W = (w/np.sqrt(1-w**2))
                model = np.linalg.lstsq(matrix * W[:,np.newaxis], x * W, rcond=None)
        except Exception as e:
            # typically, this error handled:
            # LinAlgError: SVD did not converge in Linear Least Squares
            print ('SBAS.lstsq notice:', str(e))
            return np.nan * np.zeros(matrix.shape[1])
        #print ('model', model)
        return model[0]

    def lstsq_matrix(self, pairs):
        import numpy as np
        import pandas as pd

        # also define image capture dates from interferogram date pairs 
        pairs, dates = self.pairs(pairs, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values
        
        # here are one row for every interferogram and one column for every date
        matrix = []
        for pair in pairs:
            mrow = [date>pair[0] and date<=pair[1] for date in dates]
            matrix.append(mrow)
        matrix = np.stack(matrix).astype(int)
        return matrix

    def sbas_parallel(self, pairs=None, mask=None, data='detrend', corr='corr', weight='corr',
                      chunks=None, chunksize=None, n_jobs=-1, interactive=False):

        print ('NOTE: sbas_parallel() is alias for [Weighted] Least Squares lstsq_parallel() function')

        if corr != 'corr':
            print ('NOTE: use "weight" argument instead of "corr"')

        if mask is not None:
            print ('NOTE: "mask" argument support is removed. Use "data" and "weight" arguments to call with custom data arrays')

        if chunks is not None:
            print ('NOTE: use "chunksize" argument instead of "chunks"')
            
        return self.lstsq_parallel(pairs=pairs, data=data, weight=weight,
                                   chunksize=chunksize, n_jobs=n_jobs, interactive=interactive)
                   
    def lstsq_parallel(self, pairs=None, data='detrend', weight='corr',
                       chunksize=None, n_jobs=-1, interactive=False):
        import xarray as xr
        import numpy as np
        import pandas as pd
        import dask
        from tqdm.auto import tqdm
        import joblib
        import os

        if chunksize is None:
            # smaller chunks are better for large 3D grids processing
            chunksize = self.chunksize
        #print ('chunksize', chunksize)

        # also define image capture dates from interferogram date pairs 
        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs, dates = self.pairs(pairs, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values
        
        # split large stacks to chunks about chunksize ^ 3
        minichunksize = chunksize**2 / len(pairs)
        # round to 2**deg (32,  64, 128,..., 2048)
        chunkdegrees = np.array([2**deg for deg in range(5,12) if 2**deg <= chunksize])
        minichunksize = int(chunkdegrees[np.abs(chunkdegrees-minichunksize).argmin()])
        #print ('minichunksize', minichunksize)

        # source grids lazy loading
        if isinstance(data, str):
            data = self.open_grids(pairs, data, chunksize=minichunksize, interactive=True)

        if isinstance(weight, str):
            weight = self.open_grids(pairs, weight, chunksize=minichunksize, interactive=True)
        elif isinstance(weight, (pd.Series, np.ndarray)):
            # vector weight like to average correlation
            weight = xr.DataArray(weight, dims=['pair'], coords={'pair': data.pair})

        # xarray wrapper
        model = xr.apply_ufunc(
            self.lstsq,
            data.chunk(dict(pair=-1)),
            weight.chunk(dict(pair=-1)) if weight is not None else None,
            input_core_dims=[['pair'], ['pair'] if weight is not None else []],
            exclude_dims=set(('pair',)),
            dask='parallelized',
            vectorize=True,
            output_dtypes=[np.float32],
            output_core_dims=[['date']],
            dask_gufunc_kwargs={'output_sizes': {'date': len(dates)}},
            kwargs={'matrix': self.lstsq_matrix(pairs)}
        ).rename('displacement')
        # define dates axis
        model['date'] = dates
        # set the stack index to be first
        model = model.transpose('date',...)

        if interactive:
            return model

        # use this trick to save large model using limited RAM and all CPU cores
        ts, ys, xs = model.data.blocks.shape
        assert ts == 1, 'Date chunks count should be equal to 1'
        tchunks, ychunks, xchunks = model.chunks
        coordt = model.date
        coordy = np.array_split(model.y, np.cumsum(ychunks))
        coordx = np.array_split(model.x, np.cumsum(xchunks))

        def func(iy, ix):
            chunk_filename = self.get_filenames(None, None, f'disp_chunk_{iy}_{ix}')
            if os.path.exists(chunk_filename):
                os.remove(chunk_filename)
            # lazy dask array
            data = model.data.blocks[0,iy,ix]
            # wrap the array
            da = xr.DataArray(data,
                              dims=['date','y','x'],
                              coords={'date': coordt, 'y':coordy[iy], 'x':coordx[ix]})\
                 .rename('displacement')
            # compress 3d output following the processing blocks
            netcdf_compression = self.compression()
            # use model chunks as is for better performance
            netcdf_compression['chunksizes'] = (1, minichunksize, minichunksize)
            # compute and save to NetCDF using chunks of original coordinates
            da.to_netcdf(chunk_filename,
                         unlimited_dims=['y','x'],
                         encoding={'displacement': netcdf_compression},
                         engine=self.engine,
                         compute=True)
            return chunk_filename

        # process the chunks as separate tasks on all available CPU cores
        with self.tqdm_joblib(tqdm(desc='[Correlation-Weighted] Least Squares Computing', total=ys*xs)) as progress_bar:
            filenames = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(func)(iy, ix) \
                                                     for iy in range(ys) for ix in range(xs))

        # rebuild the datasets to user-friendly format
        das = [xr.open_dataarray(f, engine=self.engine, chunks=chunksize) for f in filenames]
        das = xr.combine_by_coords(das)

        def output(dt):
            filename = self.get_filenames(None, None, f'disp_{dt}'.replace('-',''))
            if os.path.exists(filename):
                os.remove(filename)
            # compress 3d output following the processing blocks
            netcdf_compression = self.compression()
            netcdf_compression['chunksizes'] = (chunksize, chunksize)
            das.sel(date=dt).to_netcdf(filename,
                        encoding={'displacement': netcdf_compression},
                        engine=self.engine)

        # saving all the grids
        with self.tqdm_joblib(tqdm(desc='Saving', total=len(das.date))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(output)(dt) for dt in das.date.values)

        # cleanup
        for filename in filenames:
            if os.path.exists(filename):
                os.remove(filename)

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
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(ondemand)(date, dt) for (date, dt) in datetimes.items())
    
        # after merging use unmerged subswath PRM files
        # calc_dop_orb() required for SAT_baseline
        master_dt = datetimes[self.master]
        prm_ref = PRM().from_file(get_filename(master_dt)).calc_dop_orb(inplace=True)
        data = []
        for (date, dt) in datetimes.items():
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
