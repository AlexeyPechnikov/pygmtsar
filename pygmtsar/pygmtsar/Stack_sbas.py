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

class Stack_sbas(Stack_detrend):

    @staticmethod
    def lstsq1d(x, w, matrix):
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
            #print ('Stack.lstsq notice:', str(e))
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

    def lstsq(self, data, weight=None, debug=False):
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
            print ('Note: auto tune data chunk size to a half of NetCDF chunk')
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
        matrix = self.lstsq_matrix(pairs)

        def lstq_block(ys, xs):
            # 3D array
            data_block = data.isel(y=ys, x=xs).compute(n_workers=1).values.transpose(1,2,0)
            # weights can be defined by multiple ways or be set to None
            if isinstance(weight, xr.DataArray):
                # 3D array
                weight_block = weight.isel(y=ys, x=xs).compute(n_workers=1).values.transpose(1,2,0)
                # weight=1 is not allowed for the used weighted least squares calculation function 
                weight_block = np.where(weight_block>=1, 0.999999, weight_block)
                # Vectorize vec_lstsq
                vec_lstsq = np.vectorize(lambda x, w: self.lstsq1d(x, w, matrix), signature='(n),(n)->(m)')
                # Apply vec_lstsq to data_block and weight_block and revert the original dimensions order
                block = vec_lstsq(data_block, weight_block).transpose(2,0,1)
                del weight_block, vec_lstsq
            else:
                # Vectorize vec_lstsq
                vec_lstsq = np.vectorize(lambda x: self.lstsq1d(x, weight, matrix), signature='(n)->(m)')
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

    def baseline_pairs(self, days=100, meters=None, limit=None, invert=False):
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

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the sorted list of baseline pairs with reference and repeat dates,
            timelines, and baselines.

        """
        import numpy as np
        import pandas as pd

        def baseline_table():
            """
            Generates a baseline table for Sentinel-1 data, containing dates, times, and baseline components.

            This function creates a baseline table for Sentinel-1 data by processing the PRM files, which
            contain metadata for each image. The table includes dates, times, and parallel and perpendicular
            baseline components for each image.
            """
            import pandas as pd
            import numpy as np

            datetimes = self.df.datetime

            prm_ref = self.PRM()
            data = []
            for (date, dt) in datetimes.items():
                prm_rep = self.PRM(date)
                BPL, BPR = prm_ref.SAT_baseline(prm_rep).get('B_parallel', 'B_perpendicular')
                data.append({'date':date, 'BPL':BPL, 'BPR':BPR})
            df = pd.DataFrame(data).set_index('date')
            df.index = pd.DatetimeIndex(df.index)
            return df

        tbl = baseline_table()
        data = []
        for line1 in tbl.itertuples():
            counter = 0
            for line2 in tbl.itertuples():
                #print (line1, line2)
                if limit is not None and counter >= limit:
                    continue
                if not (line1.Index < line2.Index and (line2.Index - line1.Index).days < days):
                    continue
                if meters is not None and not (abs(line1.BPR - line2.BPR)< meters):
                    continue

                counter += 1
                if not invert:
                    data.append({'ref_date':line1.Index, 'rep_date': line2.Index,
                                 'ref_baseline': np.round(line1.BPR, 2),
                                 'rep_baseline': np.round(line2.BPR, 2)})
                else:
                    data.append({'ref_date':line2.Index, 'rep_date': line1.Index,
                                 'ref_baseline': np.round(line2.BPR, 2),
                                 'rep_baseline': np.round(line1.BPR, 2)})

        return pd.DataFrame(data).sort_values(['ref_date', 'rep_date'])
