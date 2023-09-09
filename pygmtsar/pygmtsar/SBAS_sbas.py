# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_detrend import SBAS_detrend
from .PRM import PRM

class SBAS_sbas(SBAS_detrend):

    @staticmethod
    def lstsq(x, w, matrix):
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
            # allow not weighted Least Squares
            w = 0.5 * np.ones_like(x)
        elif np.any(np.isnan(w)):
            # ignore pixels where correlation is not defined
            return np.nan * np.zeros(matrix.shape[1])

        assert x.shape == w.shape, f'Arrays x and w need to have equal shape, x.shape={x.shape}, w.shape={w.shape}'

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

    def stack_lstsq(self, data=None, weight=None, chunksize=None, interactive=False, debug=False):
        """
        Perform least squares (weighted or unweighted) computation on the input data in parallel.

        This function applies the least squares computation on the input data, taking into account
        interferogram date pairs, and processes it in parallel using Dask and Joblib.

        Parameters
        ----------
        pairs : array-like or pandas.DataFrame, optional
            Interferogram date pairs.
        data : str or xarray.DataArray, optional
            Input data to compute least squares on.
        weight : str, xarray.DataArray, pd.Series, or np.ndarray, optional
            Weights for the least squares computation.
        chunksize : int, optional
            Size of the chunks for processing the data.
        interactive : bool, optional
            If True, returns the result as an xarray DataArray (default is False).

        Returns
        -------
        xarray.DataArray
            The computed least squares model as an xarray DataArray.

        Examples:
        -----
        sbas.lstsq_parallel(unwraps_detrend, interactive=False)
        sbas.lstsq_parallel(unwraps_detrend, corrs, interactive=False)
        sbas.lstsq_parallel([unwraps_detrend, corrs], interactive=False)
        sbas.lstsq_parallel((unwraps_detrend, corrs), interactive=False)
        sbas.lstsq_parallel((unwraps_detrend, corrs.mean(['y', 'x'])), interactive=False)

        Notes
        -----
        This function processes large stacks by splitting them into chunks, performing the computation,
        and then rebuilding and saving the results in a user-friendly format.
    
        """
        import xarray as xr
        import numpy as np
        import pandas as pd
        import dask

        if chunksize is None:
            # 3D grids processing chunks
            # note: pair/date coordinate is not chunked
            # for base chunksize=1000 64 is slow with garbage collector notices,
            # 128 ok, 256 is faster, 512 is a bit slowly and requires more RAM
            # to have good chunks for base chunksize=512 and 1024 divider 4 looks the best   
            chunksize = self.chunksize // 4
        if debug:
            print ('DEBUG: chunksize', chunksize)

        # source grids lazy loading
        #if isinstance(data, str):
        #    data = self.open_grids(pairs, data, interactive=True)
        if data is None:
            raise ValueError('Argument data should be 3D Xarray object')
        elif weight is not None and isinstance(data, (list, tuple)):
            raise ValueErrorlueError('Weight defined twice using weight argument and as the second part in data argument')
        elif weight is None and isinstance(data, (list, tuple)):
            assert len(data) == 2, 'Argument data can be a list or a tuple of data and weight arguments'
            # this case should be processed inside lstq_block function
            data, weight   = data
            assert len(data.dims) == 3, 'Argument data should be 3D Xarray object'
        elif isinstance(data, xr.DataArray) and len(data.dims) == 3:
            # this case should be processed inside lstq_block function
            pass
        if debug:
            print ('DEBUG: data', data)

        if weight is None:
            # this case should be processed inside lstq_block function
            pass
        elif isinstance(weight, np.ndarray):
            # this case should be processed inside lstq_block function
            pass
        #elif isinstance(weight, str):
        #    weight = self.open_grids(pairs, weight, interactive=True)
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
        else:
            raise ValueError(f"Argument weight can be 1D or 3D Xarray object or Pandas Series or Numpy array or Python list")
        if debug:
            print ('DEBUG: weight', weight)

        # also define image capture dates from interferogram date pairs 
        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs, dates = self.get_pairs(data, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values
        # define pairs and dates matrix for LSQ calculation
        matrix = self.lstsq_matrix(pairs)

        def lstq_block(ys, xs):
            # 3D array
            data_block = data.isel(y=ys, x=xs).chunk(-1).compute(n_workers=1).values.transpose(1,2,0)
            # weights can be defined by multiple ways or be set to None
            if isinstance(weight, xr.DataArray):
                # 3D array
                weight_block = weight.isel(y=ys, x=xs).chunk(-1).compute(n_workers=1).values.transpose(1,2,0)
                # Vectorize vec_lstsq
                vec_lstsq = np.vectorize(lambda x, w: self.lstsq(x, w, matrix), signature='(n),(n)->(m)')
                # Apply vec_lstsq to data_block and weight_block and revert the original dimensions order
                block = vec_lstsq(data_block, weight_block).transpose(2,0,1)
                del weight_block, vec_lstsq
            else:
                # Vectorize vec_lstsq
                vec_lstsq = np.vectorize(lambda x: self.lstsq(x, weight, matrix), signature='(n)->(m)')
                # Apply vec_lstsq to data_block and weight_block and revert the original dimensions order
                block = vec_lstsq(data_block).transpose(2,0,1)
                del vec_lstsq
            del data_block
            return block

        # split to square chunks
        # use indices instead of the coordinate values to prevent the weird error raising occasionally:
        # "Reindexing only valid with uniquely valued Index objects"
        ys_blocks = np.array_split(np.arange(data.y.size), np.arange(0, data.y.size, chunksize)[1:])
        xs_blocks = np.array_split(np.arange(data.x.size), np.arange(0, data.x.size, chunksize)[1:])

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

        if interactive:
            return model

        self.save_model(model, caption='[Correlation-Weighted] Least Squares Computing', chunksize=chunksize)

    def baseline_table(self, n_jobs=-1, debug=False):
        """
        Generates a baseline table for Sentinel-1 data, containing dates, times, and baseline components.

        This function creates a baseline table for Sentinel-1 data by processing the PRM files, which
        contain metadata for each image. The table includes dates, times, and parallel and perpendicular
        baseline components for each image.

        Parameters
        ----------
        n_jobs : int, optional
            Number of CPU cores to use for parallel processing (default is -1, which means using all available cores).
        debug : bool, optional
            If True, print additional information during processing (default is False).

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the baseline table with date, times, and baseline components.

        Notes
        -----
        This function processes Sentinel-1 data by first generating PRM and LED files if they don't exist,
        then calculating doppler and orbital parameters for each image, and finally computing the baseline
        components.

        """
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
        reference_dt = datetimes[self.reference]
        prm_ref = PRM().from_file(get_filename(reference_dt)).calc_dop_orb(inplace=True)
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

    def baseline_pairs(self, days=100, meters=None, limit=None, invert=False, n_jobs=-1, debug=False):
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
            Maximum spatial separation between image pairs in meters (default is 150).
        limit : int, optional
            Maximum number of pairs to return per reference image (default is None, which means no limit).
        invert : bool, optional
            If True, invert the order of reference and repeat images (default is False).
        n_jobs : int, optional
            Number of CPU cores to use for parallel processing (default is -1, which means using all available cores).
        debug : bool, optional
            If True, print additional information during processing (default is False).

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the sorted list of baseline pairs with reference and repeat dates,
            timelines, and baselines.

        """
        import numpy as np
        import pandas as pd

        tbl = self.baseline_table(n_jobs=n_jobs, debug=debug)
        data = []
        for line1 in tbl.itertuples():
            counter = 0
            for line2 in tbl.itertuples():
                #print (line1, line2)
                if limit is not None and counter >= limit:
                    continue
                if not (line1.YDAY < line2.YDAY and line2.YDAY - line1.YDAY < days):
                    continue
                if meters is not None and not (abs(line1.BPR - line2.BPR)< meters):
                    continue

                counter += 1
                if not invert:
                    data.append({'ref_date':line1.Index, 'rep_date': line2.Index,
                                 'ref_timeline': np.round(line1.YDAY/365.25+2014, 2), 'ref_baseline': np.round(line1.BPR, 2),
                                 'rep_timeline': np.round(line2.YDAY/365.25+2014, 2), 'rep_baseline': np.round(line2.BPR, 2)})
                else:
                    data.append({'ref_date':line2.Index, 'rep_date': line1.Index,
                                 'ref_timeline': np.round(line2.YDAY/365.25+2014, 2), 'ref_baseline': np.round(line2.BPR, 2),
                                 'rep_timeline': np.round(line1.YDAY/365.25+2014, 2), 'rep_baseline': np.round(line1.BPR, 2)})

        return pd.DataFrame(data).sort_values(['ref_date', 'rep_date'])
