# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_detrend import Stack_detrend
from .PRM import PRM
import numpy as np

class Stack_sbas(Stack_detrend):

    @staticmethod
    def lstsq1d(x, w, matrix, cumsum=True):
        """
        Compute the least squares solution (or weighted least squares if weights are provided) for a given matrix.

        This function handles single-pixel processing and allows for non-weighted least squares, weighted least squares,
        and ignoring pixels where the correlation is not defined.

        Parameters
        ----------
        x : numpy.ndarray
            Input data array.
        w : numpy.ndarray or None
            Weights array for weighted least squares. If None, non-weighted least squares is used.
        matrix : numpy.ndarray
            Input matrix for which the least squares solution is computed.

        Returns
        -------
        numpy.ndarray
            Least squares solution for the given matrix.

        Notes
        -----
        If w is None, non-weighted least squares is used. If w has NaN values, the corresponding pixels are ignored.
        In case of a LinAlgError (e.g., SVD did not converge), the function returns an array of NaN values with the
        same shape as the input matrix's second dimension.

        """
        import numpy as np

        if w is None:
            # not weighted least squares calculation
            w = np.ones_like(x)
            nanmask = np.isnan(x)
        else:
            assert x.shape == w.shape, f'Arrays x and w need to have equal shape, x.shape={x.shape}, w.shape={w.shape}'
            nanmask = np.isnan(x) | np.isnan(w)
        if np.all(nanmask):
            return np.nan * np.zeros(matrix.shape[1])
        #print ('nanmask', nanmask)
        # fill nans as zeroes and set corresponding weight to 0
        # function arguments are read-only
        x = x.copy()
        # prevent weight=1
        w = (1 - 1e-6)*w.copy()
        x[nanmask] = 0
        w[nanmask] = 0
        # check if x has enough valid values
        #print (f'{x.size}, {np.sum(nanmask)}, {matrix.shape[1]}')
        #if x.size - np.sum(nanmask) < matrix.shape[1]:
        #    return np.nan * np.zeros(matrix.shape[1])

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
            #print ('Stack.lstsq notice:', str(e))
            return np.nan * np.zeros(matrix.shape[1])
        #print ('model', model)
        if not cumsum:
            return model[0]
        # mask produced cumsum zeroes by NaNs where model[0] is the timeseries values
        return np.where(~np.isnan(model[0]), np.nancumsum(model[0], axis=0), np.nan)

    def lstsq_matrix(self, pairs):
        """
        Create a matrix for use in the least squares computation based on interferogram date pairs.

        Parameters
        ----------
        pairs : pandas.DataFrame
            DataFrame containing interferogram date pairs.

        Returns
        -------
        numpy.ndarray
            Matrix with one row for every interferogram and one column for every date.
            Each element in the matrix is an integer, with 1 representing that the date
            is between the corresponding interferogram's reference and repeat dates, and
            0 otherwise.

        Notes
        -----
        This function also calculates image capture dates from the interferogram date pairs.

        """
        import numpy as np
        import pandas as pd

        # also define image capture dates from interferogram date pairs
        pairs, dates = self.get_pairs(pairs, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values

        # here are one row for every interferogram and one column for every date
        matrix = []
        for pair in pairs:
            mrow = [date>pair[0] and date<=pair[1] for date in dates]
            matrix.append(mrow)
        matrix = np.stack(matrix).astype(int)
        return matrix

    def lstsq_matrix_edge(self, pairs):
        """
        Create an edge matrix for use in the least squares computation based on interferogram date pairs.

        Parameters
        ----------
        pairs : pandas.DataFrame
            DataFrame containing interferogram date pairs.

        Returns
        -------
        numpy.ndarray
            Matrix with one row for every interferogram and one column for every date.
            Each element in the matrix is an integer, with 1 representing the end date,
            -1 the start date and 0 otherwise.

        Notes
        -----
        This function also calculates image capture dates from the interferogram date pairs.

        """
        import numpy as np
        import pandas as pd

        # also define image capture dates from interferogram date pairs
        pairs, dates = self.get_pairs(pairs, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values

        def edge(date, date1, date2):
            if date == date1:
                return -1
            if date == date2:
                return 1
            return 0
        
        # here are one row for every interferogram and one column for every date
        matrix = []
        for pair in pairs:
            mrow = [edge(date, pair[0], pair[1]) for date in dates]
            matrix.append(mrow)
        matrix = np.stack(matrix).astype(int)
        return matrix

    def lstsq(self, data, weight=None, matrix='auto', cumsum=True, debug=False):
        """
        Perform least squares (weighted or unweighted) computation on the input phase data in parallel.

        Parameters
        ----------
        data : str or xarray.DataArray, optional
            Input data to compute least squares on.
        weight : str, xarray.DataArray, pd.Series, or np.ndarray, optional
            Weights for the least squares computation.

        Returns
        -------
        xarray.DataArray
            The computed least squares model as an xarray DataArray.

        Examples:
        -----
        stack.lstsq(unwraps_detrend)
        stack.lstsq(unwraps_detrend, corrs)
        stack.lstsq(unwraps_detrend, corrs.mean(['y', 'x']))

        Notes
        -----
        This function processes large stacks by splitting them into chunks, performing the computation,
        and then rebuilding and saving the results in a user-friendly format.

        """
        import xarray as xr
        import numpy as np
        import pandas as pd
        import dask

        if isinstance(data, xr.DataArray) and len(data.dims) == 3:
            # this case should be processed inside lstq_block function
            pass
        if debug:
            print ('DEBUG: data', data)

        chunks_z, chunks_y, chunks_x = data.chunks
        if np.max(chunks_y) > self.netcdf_chunksize or np.max(chunks_x) > self.netcdf_chunksize:
            print (f'Note: data chunk size ({np.max(chunks_y)}, {np.max(chunks_x)}) is too large for stack processing')
            chunks_y = chunks_x = self.netcdf_chunksize//2
            print (f'Note: auto tune data chunk size to a half of NetCDF chunk: ({chunks_y}, {chunks_x})')
            data = data.chunk({'y': chunks_y, 'x': chunks_x})

        if weight is None:
            # this case should be processed inside lstq_block function
            pass
        elif isinstance(weight, np.ndarray):
            # this case should be processed inside lstq_block function
            pass
        elif isinstance(weight, pd.Series):
            # vector weight like to average correlation
            weight = weight.values
        elif isinstance(weight, (list, tuple)):
            weight = np.asarray(weight)
        elif isinstance(weight, xr.DataArray) and len(weight.dims) == 1:
            # 1D array
            weight = weight.values
        elif isinstance(weight, xr.DataArray) and len(weight.dims) == 3:
            # this case should be processed inside lstq_block function
            assert weight.shape == data.shape, 'ERROR: data and weight dataarrays should have the same dimensions'
            if np.max(chunks_y) > self.netcdf_chunksize or np.max(chunks_x) > self.netcdf_chunksize:
                print ('Note: auto tune weight chunk size to a half of NetCDF chunk')
                weight = weight.chunk({'y': chunks_y, 'x': chunks_x})
        else:
            raise ValueError(f"Argument weight can be 1D or 3D Xarray object or Pandas Series or Numpy array or Python list")
        if debug:
            print ('DEBUG: weight', weight)

        # also define image capture dates from interferogram date pairs 
        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs, dates = self.get_pairs(data, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values
        
        # define pairs and dates matrix for LSQ calculation
        if isinstance(matrix, str) and matrix == 'auto':
            if debug:
                print ('DEBUG: Generate default matrix')
            matrix = self.lstsq_matrix(pairs)
        else:
            if debug:
                print ('DEBUG: Use user-supplied matrix')

        def lstq_block(ys, xs):
            # 3D array
            data_block = data.isel(y=ys, x=xs).compute(n_workers=1).values.transpose(1,2,0)
            # weights can be defined by multiple ways or be set to None
            if isinstance(weight, xr.DataArray):
                # 3D array
                weight_block = weight.isel(y=ys, x=xs).compute(n_workers=1).values.transpose(1,2,0)
                # weight=1 is not allowed for the used weighted least squares calculation function 
                weight_block = np.where(weight_block>=1, 1, weight_block)
                # Vectorize vec_lstsq
                vec_lstsq = np.vectorize(lambda x, w: self.lstsq1d(x, w, matrix, cumsum), signature='(n),(n)->(m)')
                # Apply vec_lstsq to data_block and weight_block and revert the original dimensions order
                block = vec_lstsq(data_block, weight_block).transpose(2,0,1)
                del weight_block, vec_lstsq
            else:
                # Vectorize vec_lstsq
                vec_lstsq = np.vectorize(lambda x: self.lstsq1d(x, weight, matrix, cumsum), signature='(n)->(m)')
                # Apply vec_lstsq to data_block and weight_block and revert the original dimensions order
                block = vec_lstsq(data_block).transpose(2,0,1)
                del vec_lstsq
            del data_block
            return block

        # split to chunks
        # use indices instead of the coordinate values to prevent the weird error raising occasionally:
        # "Reindexing only valid with uniquely valued Index objects"
        # re-check the chunk sizes as it can be tunned above
        chunks_z, chunks_y, chunks_x = data.chunks
        ys_blocks = np.array_split(np.arange(data.y.size), np.cumsum(chunks_y)[:-1])
        xs_blocks = np.array_split(np.arange(data.x.size), np.cumsum(chunks_x)[:-1])

        blocks_total = []
        for ys_block in ys_blocks:
            blocks = []
            for xs_block in xs_blocks:
                block = dask.array.from_delayed(dask.delayed(lstq_block)(ys_block, xs_block),
                                                shape=(len(dates), ys_block.size, xs_block.size),
                                                dtype=np.float32)
                blocks.append(block)
                del block
            blocks_total.append(blocks)
            del blocks
        model = dask.array.block(blocks_total)
        del blocks_total
        coords = {'date': pd.to_datetime(dates), 'y': data.y.values, 'x': data.x.values}
        model = xr.DataArray(model, coords=coords).rename('displacement')

        return model
        #self.save_cube(model, caption='[Correlation-Weighted] Least Squares Computing')

    def baseline_pairs(self, days=100, meters=None, invert=False, **kwargs):
        print('Note: function baseline_pairs() renamed to sbas_pairs(). Use separate filtering functions when needed.')
        return self.sbas_pairs(days=days, meters=meters, invert=invert)
    
    def sbas_pairs_filter_dates(self, pairs, dates):
        return pairs[(~pairs['ref'].isin(dates))&(~pairs['rep'].isin(dates))]

    def sbas_pairs_limit(self, pairs, limit=2, iterations=1):
        """
        min_limit : int, optional
                Minimum number of pairs per date (default is None, which means no limit).
        max_limit : int, optional
            Maximum number of pairs per date (default is None, which means no limit).    
        """
        import numpy as np

        df = pairs.copy()

        # extend simplified definition
        if isinstance(limit, int):
            limit = (limit, None)

        if limit[0] is not None:
            # clean hanging nodes
            for lim in np.repeat(range(1, limit[0] + 1), iterations):
                #print ('lim', lim)
                dates, counts = np.unique(df.values[:,:2].reshape(-1), return_counts=True)
                dates = dates[counts>=lim]
                df = df[(df['ref'].isin(dates))&(df['rep'].isin(dates))]

        if limit[1] is not None:
            print ('TODO: process upper limit')

        return df

    def sbas_pairs(self, days=100, meters=None, invert=False, dates=None):
        """
        Generates a sorted list of baseline pairs based on specified temporal and spatial criteria.
    
        This function creates a list of baseline pairs for Sentinel-1 data, considering the given
        temporal and spatial constraints. The list includes pairs that meet the specified days and
        meters criteria.
    
        Parameters
        ----------
        days : int, optional
            Maximum temporal separation between image pairs in days (default is 100).
        meters : int, optional
            Maximum spatial separation between image pairs in meters (default is None).
        invert : bool, optional
            If True, invert the order of reference and repeat images (default is False).
    
        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the sorted list of baseline pairs with reference and repeat dates,
            timelines, and baselines.
    
        """
        import numpy as np
        import pandas as pd
    
        def baseline_table(dates):
            """
            Generates a baseline table for Sentinel-1 data, containing dates, times, and baseline components.
    
            This function creates a baseline table for Sentinel-1 data by processing the PRM files, which
            contain metadata for each image. The table includes dates, times, and parallel and perpendicular
            baseline components for each image.
            """
            import pandas as pd
            import numpy as np
    
            prm_ref = self.PRM()
            data = []
            for date in dates:
                prm_rep = self.PRM(date)
                BPL, BPR = prm_ref.SAT_baseline(prm_rep).get('B_parallel', 'B_perpendicular')
                data.append({'date':date, 'BPL':BPL, 'BPR':BPR})
            df = pd.DataFrame(data).set_index('date')
            df.index = pd.DatetimeIndex(df.index)
            return df
    
        if dates is None:
            dates = self.df.index
    
        tbl = baseline_table(dates)
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

    def sbas_pairs_extend(self, baseline_pairs):
        import pandas as pd
        import numpy as np

        refreps = [(key[0], key[1]) for key in baseline_pairs[['ref', 'rep']].values]
        pairs = []
        for _, row1 in baseline_pairs.iterrows():
            for _, row2 in baseline_pairs.iterrows():
                if row1['ref'] >= row2['rep']:
                    continue
                pair = (row1['ref'], row2['rep'], row1['rep'], row1['ref_baseline'], row2['rep_baseline'])
                key = (row1['ref'], row2['rep'])
                if key in refreps:
                    continue
                refreps.append(key)
                pairs.append(pair)
        df = pd.DataFrame(pairs, columns=['ref', 'rep', 'rel', 'ref_baseline', 'rep_baseline'])
        return pd.concat([baseline_pairs.assign(rel=np.datetime64('nat')),
                          df.assign(pair=[f'{ref} {rep}' for ref, rep in zip(df['ref'].dt.date, df['rep'].dt.date)],
                                    baseline=df.rep_baseline - df.ref_baseline,
                                    duration=(df['rep'] - df['ref']).dt.days)]).sort_values(['ref', 'rep'])

    @staticmethod
    def sbas_pairs_covering(pairs, column, count, func='min'):
        import pandas as pd

        df = pairs.copy()
        # Generate the new "date" column as a list of dates between 'ref' and 'rep' dates
        df['date'] = df.apply(lambda row: pd.date_range(row['ref'], row['rep']).tolist(), axis=1)

        # Expand the list into individual rows in the new dataframe
        df = df.explode('date').reset_index(drop=True)

        df_grouped = df.groupby('date')[column]
        # filter by lowest/largest per date
        if func == 'min':
            df_selected = df_grouped.nsmallest(count)
        elif func == 'max':
            df_selected = df_grouped.nlargest(count)
        else:
            raise ValueError(f"Unsupported function {func}. Should be 'min' or'max'")
        df = df.merge(df_selected.reset_index(level=1)['level_1'].reset_index(), left_index=True, right_on='level_1')
        df = df.drop(columns=['level_1']).rename(columns={'date_x': 'date'})
        del df['date_y'], df['date']
        # drop duplicates
        df = df.drop_duplicates(subset=['ref', 'rep']).reset_index(drop=True)

        return df

    def sbas_pairs_covering_deviation(self, pairs, count, column='stddev'):
        return self.sbas_pairs_covering(pairs, column, count, 'min')

    def sbas_pairs_covering_correlation(self, pairs, count, column='corr'):
        return self.sbas_pairs_covering(pairs, column, count, 'max')

    def sbas_pairs_extend(self, baseline_pairs):
        import pandas as pd
        import numpy as np

        refreps = [(key[0], key[1]) for key in baseline_pairs[['ref', 'rep']].values]
        pairs = []
        for _, row1 in baseline_pairs.iterrows():
            for _, row2 in baseline_pairs.iterrows():
                if row1['ref'] >= row2['rep']:
                    continue
                if row1['rep'] != row2['ref']:
                    continue
                pair = (row1['ref'], row2['rep'], row1['rep'], row1['ref_baseline'], row2['rep_baseline'])
                key = (row1['ref'], row2['rep'])
                if key in refreps:
                    continue
                refreps.append(key)
                pairs.append(pair)
        df = pd.DataFrame(pairs, columns=['ref', 'rep', 'rel', 'ref_baseline', 'rep_baseline'])
        return pd.concat([baseline_pairs.assign(rel=np.datetime64('nat')),
                          df.assign(pair=[f'{ref} {rep}' for ref, rep in zip(df['ref'].dt.date, df['rep'].dt.date)],
                                    baseline=df.rep_baseline - df.ref_baseline,
                                    duration=(df['rep'] - df['ref']).dt.days)]).sort_values(['ref', 'rep'])

    def correlation_extend(self, corr, baseline_pairs):
        import xarray as xr
        import pandas as pd
    
        corr_rels = []
        for (ref, rel, rep) in zip(baseline_pairs.ref,
                                   baseline_pairs.rel,
                                   baseline_pairs.rep):
            #print (ref, rel, rep)
            if pd.isnull(rel):
                #print (corr.sel(pair=f'{ref} {rep}'))
                corr_rels.append(corr.sel(pair=f'{ref.date()} {rep.date()}'))
            else:
                corr_ref = corr.sel(pair=[f'{ref.date()} {rel.date()}'])
                corr_rep = corr.sel(pair=[f'{rel.date()} {rep.date()}'])
                corr_rel = xr.concat([corr_ref, corr_rep], dim='pair').min('pair')\
                    .assign_coords(pair=f'{ref.date()} {rep.date()}', ref=ref, rep=rep)
                #print (corr_rel)
                corr_rels.append(corr_rel)
                del corr_ref, corr_rep, corr_rel
        return xr.concat(corr_rels, dim='pair')

    def interferogram_extend(self, phase, baseline_pairs, wrap=True):
        import xarray as xr
        import pandas as pd
    
        phase_rels = []
        for (ref, rel, rep) in zip(baseline_pairs.ref,
                                   baseline_pairs.rel,
                                   baseline_pairs.rep):
            #print (ref, rel, rep)
            if pd.isnull(rel):
                #print (phase.sel(pair=f'{ref} {rep}'))
                phase_rels.append(phase.sel(pair=f'{ref.date()} {rep.date()}'))
            else:
                phase_ref = phase.sel(pair=[f'{ref.date()} {rel.date()}'])
                phase_rep = phase.sel(pair=[f'{rel.date()} {rep.date()}'])
                phase_rel = xr.concat([phase_ref, phase_rep], dim='pair').sum('pair')\
                    .assign_coords(pair=f'{ref.date()} {rep.date()}', ref=ref, rep=rep)
                #print (phase_rel)
                phase_rels.append(phase_rel)
                del phase_ref, phase_rep, phase_rel
        if wrap:
            return self.wrap(xr.concat(phase_rels, dim='pair'))
        return xr.concat(phase_rels, dim='pair')

    def sbas_pairs_fill(self, data):
        """
        sbas.sbas_pairs_fill(pairs_best)
        """
        import pandas as pd
        import numpy as np

        baseline_pairs = self.get_pairs(data)
    
        matrix = self.lstsq_matrix(baseline_pairs).max(axis=0)
        #print ('matrix', matrix)
        dates = self.get_pairs(baseline_pairs, dates=True)[1]
        #print ('dates', dates)
        # ignore the first date
        pairs = []
        for date in dates[matrix==0][1:]:
            date_idx = np.where(dates==date)[0][0]
            #print ('np.where(dates==date)', np.where(dates==date))
            date_prev = dates[date_idx-1]
            #print (date_prev, '=>', date)
            pairs.append([pd.Timestamp(date_prev), pd.Timestamp(date)])
        if pairs == []:
            return
    
        df = pd.DataFrame(pairs, columns=['ref', 'rep'])
        return df.assign(pair=[f'{ref} {rep}' for ref, rep in zip(df['ref'].dt.date, df['rep'].dt.date)],
                                    duration=(df['rep'] - df['ref']).dt.days).sort_values(['ref', 'rep'])

    def interferogram_fill(self, phase):
        """
        phase_sbas_fill = sbas.interferogram_fill(phase_sbas.sel(pair=phase_sbas.pair.isin(pairs_best)))
        """
        import xarray as xr
        import pandas as pd
    
        pairs = self.sbas_pairs_fill(phase)
        if pairs is None:
            return phase

        empty = xr.zeros_like(phase.isel(pair=0))
        phases = [phase]
        for (ref, rep) in zip(pairs.ref, pairs.rep):
            #print (ref, rep)
            phases.append(empty.assign_coords(pair=f'{ref.date()} {rep.date()}', ref=ref, rep=rep))
        return xr.concat(phases, dim='pair')

    def correlation_fill(self, corr):
        """
        corr_sbas_fill = sbas.correlation_fill(corr_sbas_ext.sel(pair=corr_sbas_ext.pair.isin(pairs_best)))
        """
        import xarray as xr
        import pandas as pd
    
        pairs = self.sbas_pairs_fill(corr)
        if pairs is None:
            return corr

        empty = xr.ones_like(corr.isel(pair=0))
        corrs = [corr]
        for (ref, rep) in zip(pairs.ref, pairs.rep):
            #print (ref, rep)
            corrs.append(empty.assign_coords(pair=f'{ref.date()} {rep.date()}', ref=ref, rep=rep))
        return xr.concat(corrs, dim='pair')

    def baseline_plot(self, baseline_pairs):
        print ('NOTE: this function is deprecated, use instead Stack.plot_baseline()')
        self.plot_baseline(baseline_pairs)
        
    def plot_baseline(self, baseline_pairs):
        import numpy as np
        import pandas as pd
        import seaborn as sns
        import adjustText
        import matplotlib.pyplot as plt

        plt.figure()

        # plot dates/baselines marks
        df = pd.DataFrame(np.concatenate([baseline_pairs[['ref', 'ref_baseline']],
                                     baseline_pairs[['rep', 'rep_baseline']]]),
                    columns=['date', 'baseline']).drop_duplicates()
        sns.scatterplot(x='date', y='baseline', data=df, marker='o', color='b', s=40)
        # plot reference date on top
        sns.scatterplot(x='date', y='baseline', data=df[df.date==self.reference], marker='o', color='r', s=40, zorder=1000)

        # plot pairs
        for _, row in baseline_pairs.iterrows():
            plt.plot([row['ref'], row['rep']], [row['ref_baseline'], row['rep_baseline']],
                     c='#30a2da' if pd.isnull(row['rel']) else '#2dda30', lw=0.5)
        # highlight the longest pair
        # for _, row in self.get_pairs(baseline_pairs).sort_values('duration', ascending=False).head(1).iterrows():
        #     plt.plot([row['ref'], row['rep']], [row['ref_baseline'], row['rep_baseline']], c='red', lw=1)

        # Create annotations with adjust_text
        texts = []
        for x, y in df.values:
            texts.append(plt.text(x, y, str(x.date()), ha='center', va='bottom',
                                  c='r' if str(x.date()) == self.reference else 'black'))
        adjustText.adjust_text(texts)

        plt.xlabel('Timeline')
        plt.ylabel('Perpendicular Baseline, [m]')
        plt.title('Baseline')
        plt.grid()

    def plot_baseline_duration(self, baseline_pairs, interval_days=6, caption='Durations Histogram',
                               column=None, ascending=None, cmap='turbo', vmin=None, vmax=None):
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        max_duration = baseline_pairs.duration.max()
        bins = np.arange(interval_days / 2, max_duration + interval_days, interval_days)
        bin_midpoints = (bins[:-1] + bins[1:]) / 2
        #print ('bins', len(bins), bins)

        fig, ax = plt.subplots()

        if column is not None and ascending is None:
            # Calculate histogram with average column values
            counts, edges = np.histogram(baseline_pairs.duration, bins=bins)
            sums, _ = np.histogram(baseline_pairs.duration, bins=bins, weights=baseline_pairs[column])
            averages = sums / counts

            # Normalize the average values for coloring
            norm = mcolors.Normalize(vmin=vmin if vmin is not None else np.nanmin(averages),
                                     vmax=vmax if vmax is not None else np.nanmax(averages))
            cmap = plt.cm.get_cmap(cmap)

            for i in range(len(bin_midpoints)):
                bin_color = 'white' if np.isnan(averages[i]) else cmap(norm(averages[i]))
                ax.bar(bin_midpoints[i], counts[i], width=bins[i+1] - bins[i], color=bin_color,
                       edgecolor='black', align='center', zorder=3)

            plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label=f'Average {column}')
        elif column is not None:
            norm = mcolors.Normalize(vmin=vmin if vmin is not None else baseline_pairs[column].min(),
                                     vmax=vmax if vmax is not None else baseline_pairs[column].max())
            cmap = plt.cm.get_cmap(cmap)

            for i in range(len(bins) - 1):
                bin_data = baseline_pairs[(baseline_pairs.duration >= bins[i]) & (baseline_pairs.duration < bins[i + 1])]
                bin_data = bin_data.sort_values(by=column, ascending=ascending)
                #print (i, bin_data[column].mean())

                bottom = 0
                for _, row in bin_data.iterrows():
                    #print (i, row[column])
                    color = cmap(norm(row[column]))
                    ax.bar(bin_midpoints[i], 1, bottom=bottom, width=bins[i+1] - bins[i],
                           color=color, edgecolor='black', align='center', zorder=3)
                    bottom += 1

            plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label=f'{column} value')
        else:
            ax.hist(baseline_pairs.duration, bins=bins, color='skyblue', edgecolor='black', align='mid', zorder=3)

        ax.set_xlabel('Duration [days]')
        ax.set_ylabel('Count')
        ax.set_title(caption)
        ax.set_xticks(bin_midpoints if interval_days >= 12 else bin_midpoints[1::2])
        ax.set_xlim(interval_days / 2, max_duration + interval_days / 2)
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_horizontalalignment('center')
        ax.grid(True, color='lightgrey', zorder=0)

    def plot_baseline_attribute(self, baseline_pairs, pairs_best=None, column=None, caption='Baseline Attribute'):
        import numpy as np
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        assert column is not None, 'ERROR: specify the column argument'

        plt.figure()

        # plot dates/baselines marks
        df = pd.DataFrame(np.concatenate([baseline_pairs[['ref', column]], baseline_pairs[['rep', column]]]),
                                         columns=['date', column]).drop_duplicates()
        sns.scatterplot(x='date', y=column, data=df, marker='o', color='r', s=10, label='All Pairs')

        if pairs_best is not None:
            df_best = pd.DataFrame(np.concatenate([pairs_best[['ref', column]], pairs_best[['rep', column]]]),
                                             columns=['date', column]).drop_duplicates()
            sns.scatterplot(x='date', y=column, data=df_best, marker='o', color='g', s=40, label='Selected Pairs')

        # plot pairs
        for _, row in baseline_pairs.iterrows():
            plt.plot([row['ref'], row['rep']], [row[column], row[column]], c='r', lw=0.5)

        if pairs_best is not None:
            for _, row in pairs_best.iterrows():
                plt.plot([row['ref'], row['rep']], [row[column], row[column]], c='g', lw=1)

        plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncols=2)
        plt.xlabel('Timeline')
        plt.ylabel(f'Column "{column}"')
        plt.title(caption, y=1.2)
        plt.grid()

    def plot_baseline_correlation(self, baseline_pairs, pairs_best=None, column='corr'):
        #print ('NOTE: this function is deprecated, use instead Stack.plot_baseline_attribute()')
        self.plot_baseline_attribute(baseline_pairs, pairs_best, column=column, caption='Baseline Correlation')

    def plot_baseline_deviation(self, baseline_pairs, pairs_best=None, column='stddev'):
        #print ('NOTE: this function is deprecated, use instead Stack.plot_baseline_attribute()')
        self.plot_baseline_attribute(baseline_pairs, pairs_best, column=column, caption='Baseline Deviation')

    def plot_baseline_displacement(self, phase, corr=None, caption=None, cmap='turbo',
                                   displacement=True, unwrap=True,
                                   stl=False, stl_freq='W', stl_periods=52, stl_robust=True,
                                   los=False, tolerance=np.pi/2, xlabel_rotation=45,
                                   legend=True, legend_alpha=None):
        """
        Performs 1D unwrapping, linear regression, and STL on a given set of phase values.
    
        The linear regression model is represented as:
            y = β0 + β1 * x
    
        Where:
            y: Dependent variable (the outcome being predicted).
            x: Independent variable (the predictor).
            β0: Intercept (the value of y when x is 0).
            β1: Slope or "velocity" (the rate of change in y for a one-unit change in x).
    
        In this model, 'β' (beta) symbols followed by indices (0, 1, 2, ...) represent 
        the coefficients in the regression equation. β0 is always the intercept, and 
        β1, β2, etc., represent the coefficients of the predictor variables.
    
        """
        import numpy as np
        import xarray as xr
        import pandas as pd
        from scipy.stats import linregress
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
    
        if not displacement and stl:
            print ("NOTE: Displacement is automatically set to 'True' because it is required for 'stl=True'.")
        assert isinstance(phase, xr.DataArray) and phase.dims == ('pair',), \
            'ERROR: Argument phase should be 1D Xarray with "pair" dimension'
        plt.figure()
        colors = matplotlib.cm.get_cmap(cmap)
    
        df = phase.to_dataframe()
        df['corr'] = corr.values if corr is not None else 1
        pairs, dates = self.get_pairs(phase, dates=True)
        dates = pd.DatetimeIndex(dates)
        matrix = self.lstsq_matrix(pairs)
    
        if unwrap:
            df['phase'] = self.unwrap_pairs(phase.values, df['corr'].values, matrix, tolerance)
        
        unit = 'rad'
        name = 'Phase'
        if los:
            df['phase'] = self.los_displacement_mm(df['phase'])
            unit = 'mm'
            name = 'Displacement'
    
        if displacement or stl:
            solution = self.lstsq1d(df['phase'].values, 0.999*df['corr'].values if corr is not None else None, matrix)
            #print ('solution', solution)
            days = (dates - dates[0]).days
            #print ('days', days)
            slope, intercept, r_value, p_value, std_err = linregress(days, solution)
            #print (slope, intercept, r_value, p_value, std_err)
            velocity = np.round(slope*365.25, 2)
            values = intercept + slope*days
    
        if stl:
            #assert displacement, 'ERROR: Argument value stl=True requires argument displacement=True'
            dt, dt_periodic = self.stl_periodic(dates, stl_freq)
            trend, seasonal, resid = self.stl1d(solution, dt, dt_periodic, periods=stl_periods, robust=stl_robust)
            #years = ((dt_periodic.date[-1] - dt_periodic.date[0]).dt.days/365.25).item()
            stl_dates = pd.to_datetime(dt_periodic)
            stl_days = (stl_dates - stl_dates[0]).days
            #print ('stl_days', stl_days)
            stl_slope, stl_intercept, stl_r_value, stl_p_value, stl_std_err = linregress(stl_days, trend)
            #print (stl_slope, stl_intercept, stl_r_value, stl_p_value, stl_std_err)
            stl_velocity = np.round(stl_slope*365.25, 2)
            # for linear trend
            #velocity = (trend[-1] - trend[0])/years
            # for integral trend
            #velocity = np.round(years*np.nansum(trend - np.nanmin(trend))/trend.size, 2)
    
        vmin = df['phase'].min()
        vmax = df['phase'].max()
        y_min = y_max = 0
        dts = []
        idx = -1
        errors = []
        for i, row in df.iterrows():
            idx += 1
            for item in ['ref', 'rep']:
                if row[item] not in dts:
                    plt.axvline(x=row[item], color='black', linewidth=0.5, linestyle='--')
                    dts.append(row[item])
            if displacement or stl:
                # find pair fererence date timeline position
                position1 = np.where(np.asarray(dates)==row['ref'])
                position2 = np.where(np.asarray(dates)==row['rep'])
                start = solution[position1][0]
                end = solution[position2][0]
                #print (row['phase'], (end - start), row['phase'] - (end - start))
                errors.append(row['phase'] - (end - start))
            else:
                start = 0
            y_min = min(y_min, start, row['phase'] + start)
            y_max = max(y_max, start, row['phase'] + start)
            plt.plot([row['ref'], row['rep']], [start, row['phase'] + start],
                     c=colors(row['corr']) if corr is not None else 'grey', alpha=0.8,
                     linewidth=0.5)
        #print ('errors', errors)
        #print ('y_min', y_min, 'y_max', y_max)
        if displacement or stl:
            errors = np.asarray(errors)
            #print ('errors', errors)
            weights = df['corr'].values
            nanmask = np.isnan(errors) | np.isnan(weights)
            errors = errors[~nanmask]
            weights = weights[~nanmask]
            #print ('weights', weights)
            #print (errors.size, weights.size, errors)
            #rmse = np.sqrt(np.sum(errors**2) / errors.size)
            rmse = np.sqrt(np.sum(weights * (errors**2)) / np.sum(weights) / errors.size)
            #print ('weighted PI-scaled rmse', np.round(rmse / np.pi, 2))
            lsq_mean = np.nanmean(solution)
            plt.plot(dates, solution, color='black', linestyle='--', linewidth=2, label=f'LSQ mean={lsq_mean:0.2f} [{unit}]')
            plt.axhline(y=lsq_mean, color='black', linestyle=':')
            plt.plot(dates, values, color='blue', linestyle='-', linewidth=2,
                     label=f'LSQ β1={velocity:0.1f} and β0={intercept:0.1f} [{unit}/year], P-value={p_value:0.2f}')
    
        if stl:
            plt.plot(dt_periodic.date, trend, color='blue', linestyle='--', linewidth=2,
                     label=f'STL β1={stl_velocity:0.1f} and β0={stl_intercept:0.1f} [{unit}/year]')
            plt.plot(dt_periodic.date, seasonal, color='green', linestyle='--', linewidth=1, label='STL Seasonal')
            plt.plot(dt_periodic.date, resid, color='red', linestyle='--', linewidth=1, label='STL Residual')
    
        plt.grid(axis='y')
        if unit == 'rad':
            y_min = np.floor(y_min / np.pi) * np.pi
            y_max = np.ceil(y_max / np.pi) * np.pi
            y_ticks = np.arange(y_min, y_max + np.pi, np.pi)
            plt.yticks(y_ticks, [f'{y:.0f}π' if y != 0 else '0' for y in y_ticks / np.pi])
    
        if corr is not None:
            sm = ScalarMappable(cmap=colors, norm=plt.Normalize(vmin=0, vmax=1))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=plt.gca(), orientation='vertical')
            cbar.set_label('Correlation')
    
        plt.xticks(rotation=xlabel_rotation)
        plt.xlabel('Timeline')
        plt.ylabel(f'{name}, [{unit}]')
        plt.title('Displacement' \
                  + (f' RMSE={rmse:0.3f} [{unit}]' if displacement or stl else '') \
                  + (f'\n{caption}' if caption is not None else ''))
        plt.xlim([dates[0], dates[-1]])
        if (displacement or stl) and legend:
            plt.legend(framealpha=legend_alpha)

    def plot_baseline_displacement_los_mm(self, phase, corr=None, caption=None, cmap='turbo',
                                   displacement=True, unwrap=True,
                                   stl=False, stl_freq='W', stl_periods=52, stl_robust=True,
                                   tolerance=np.pi/2, xlabel_rotation=45,
                                   legend=True, legend_alpha=None):
        self.plot_baseline_displacement(phase=phase, corr=corr, caption=caption, cmap=cmap,
                                   displacement=displacement, unwrap=unwrap,
                                   stl=stl, stl_freq=stl_freq, stl_periods=stl_periods, stl_robust=stl_robust,
                                   los=True, tolerance=tolerance, xlabel_rotation=xlabel_rotation,
                                   legend=legend, legend_alpha=legend_alpha)

    def rmse(self, data, solution, weight=None):
        """
        Calculate difference between pairs and dates
    
        Use to calculate solution vs pair unwrapped phases difference as
        diff = sbas.stack_vs_cube(phase_unwrap, solution)
        error = np.sqrt(sbas.wrap(diff)**2).sum('pair')
        """
        import numpy as np
        import xarray as xr
    
        # extract pairs
        pairs = self.get_pairs(data)
        # unify data and solution
        pairs = pairs[pairs.ref.isin(solution.date.values)&pairs.rep.isin(solution.date.values)]
        # calculate differences between end and start dates for all the pairs
        error_pairs = []
        for rec in pairs.itertuples():
            #print (rec.rep, rec.ref, rec.Index, rec)
            error_pair = self.wrap(data.sel(pair=rec.pair) - (solution.sel(date=rec.rep) - solution.sel(date=rec.ref)))
            error_pairs.append(error_pair**2)
        # form 3D stack
        error = xr.concat(error_pairs, dim='pair').assign_coords({'pair': pairs.pair})
        return np.sqrt((weight * error).sum('pair') / weight.sum('pair') / len(pairs))

    def plot_displacement(self, data, caption='Cumulative LOS Displacement, [rad]',
                          quantile=None, vmin=None, vmax=None, symmetrical=False, aspect=None, **kwargs):
        import numpy as np
        import matplotlib.pyplot as plt

        if quantile is not None:
            assert vmin is None and vmax is None, "ERROR: arguments 'quantile' and 'vmin', 'vmax' cannot be used together"

        if quantile is not None:
            vmin, vmax = np.nanquantile(data, quantile)

        # define symmetrical boundaries
        if symmetrical is True and vmax > 0:
            minmax = max(abs(vmin), vmax)
            vmin = -minmax
            vmax =  minmax

        plt.figure()
        data.plot.imshow(vmin=vmin, vmax=vmax, cmap='turbo')
        self.plot_AOI(**kwargs)
        self.plot_POI(**kwargs)
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        plt.title(caption)

    def plot_displacements(self, data, caption='Cumulative LOS Displacement, [rad]', cols=4, size=4, nbins=5, aspect=1.2, y=1.05,
                           quantile=None, vmin=None, vmax=None, symmetrical=False):
        import numpy as np
        import matplotlib.pyplot as plt

        if quantile is not None:
            assert vmin is None and vmax is None, "ERROR: arguments 'quantile' and 'vmin', 'vmax' cannot be used together"

        if quantile is not None:
            vmin, vmax = np.nanquantile(data, quantile)

        # define symmetrical boundaries
        if symmetrical is True and vmax > 0:
            minmax = max(abs(vmin), vmax)
            vmin = -minmax
            vmax =  minmax

        # multi-plots ineffective for linked lazy data
        fg = data.plot.imshow(
            col='date',
            col_wrap=cols, size=size, aspect=aspect,
            vmin=vmin, vmax=vmax, cmap='turbo'
        )
        fg.set_axis_labels('Range', 'Azimuth')
        fg.set_ticks(max_xticks=nbins, max_yticks=nbins)
        fg.fig.suptitle(caption, y=y)

    def plot_velocity(self, data, caption='Velocity, mm/year',
                      quantile=None, vmin=None, vmax=None, symmetrical=False, aspect=None, alpha=1, **kwargs):
        import numpy as np
        import matplotlib.pyplot as plt
    
        if quantile is not None:
            assert vmin is None and vmax is None, "ERROR: arguments 'quantile' and 'vmin', 'vmax' cannot be used together"
    
        if quantile is not None:
            vmin, vmax = np.nanquantile(data, quantile)
    
        # define symmetrical boundaries
        if symmetrical is True and vmax > 0:
            minmax = max(abs(vmin), vmax)
            vmin = -minmax
            vmax =  minmax
    
        plt.figure()
        data.plot.imshow(vmin=vmin, vmax=vmax, alpha=alpha, cmap='turbo')
        self.plot_AOI(**kwargs)
        self.plot_POI(**kwargs)
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        plt.title(caption)

    def plot_rmse(self, rmse, caption='RMSE', cmap='turbo',
                  quantile=None, vmin=None, vmax=None, symmetrical=False, **kwargs):
        import numpy as np
        import matplotlib.pyplot as plt
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')

        if quantile is not None:
            assert vmin is None and vmax is None, "ERROR: arguments 'quantile' and 'vmin', 'vmax' cannot be used together"

        if quantile is not None:
            vmin, vmax = np.nanquantile(rmse, quantile)

        # define symmetrical boundaries
        if symmetrical is True and vmax > 0:
            minmax = max(abs(vmin), vmax)
            vmin = -minmax
            vmax =  minmax

        plt.figure()
        rmse.plot.imshow(cmap=cmap, vmin=vmin, vmax=vmax)
        self.plot_AOI(**kwargs)
        self.plot_POI(**kwargs)
        plt.title(caption)
