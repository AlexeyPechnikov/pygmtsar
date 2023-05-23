# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_unwrap import SBAS_unwrap

class SBAS_detrend(SBAS_unwrap):

    def detrend(self, dataarray, fit=None, fit_intercept=True, fit_dem=True, fit_coords=True,
            resolution_meters=90, wavelength=None, truncate=3.0, debug=False):
        """
        Detrend the input data array by combining optional topography and linear components removal and Gaussian filtering.

        Parameters
        ----------
        dataarray : xarray.DataArray
            The input data array.
        fit : bool, optional
            Whether to apply the same detrend options to all components (fit_intercept, fit_dem, fit_coords).
        fit_intercept : bool, optional
            Whether to remove the mean value (plane).
        fit_dem : bool, optional
            Whether to detrend the topography.
        fit_coords : bool, optional
            Whether to detrend the linear coordinate components.
        resolution_meters : float, optional
            The processing resolution in meters to prevent overfitting and reduce grid size.
        wavelength : float, tuple, list, or ndarray, optional
            The cut-off wavelength(s) for the Gaussian filter in meters.
            If a tuple, list, or ndarray is provided, the output is obtained by subtracting the
            filtered data arrays corresponding to the minimum and maximum wavelength values.
            If None, no filtering is performed.
        truncate : float, optional
            The filter window size in sigmas.
        debug : bool, optional
            Whether to print debug information.

        Returns
        -------
        xarray.DataArray
            The detrended data array.

        Examples
        --------
        Simplest detrending:
        sbas.detrend(pair)

        Detrend the unwrapped interferogram in radar coordinates.
        See the following GitHub issues for more details:
        - https://github.com/gmtsar/gmtsar/issues/98
        - https://github.com/gmtsar/gmtsar/issues/411
        """
        import numpy as np

        if fit is not None:
            # set all fir options the same
            fit_intercept = fit_dem = fit_coords = fit

        if wavelength is None:
            out1 = dataarray
        elif isinstance(wavelength, (tuple, list, np.ndarray)):
            # range-pass
            assert len(wavelength) == 2, 'ERROR: wavelength argument should be a list of two elements or scalar on None (omitted)'
            out1 = self._gaussian(dataarray, wavelength=np.min(wavelength), truncate=truncate,
                                    resolution_meters=resolution_meters, debug=debug) \
                 - self._gaussian(dataarray, wavelength=np.max(wavelength), truncate=truncate,
                                    resolution_meters=resolution_meters, debug=debug)
        else:
            # high-pass
            out1 = dataarray - self._gaussian(dataarray, wavelength=wavelength, truncate=truncate,
                     resolution_meters=resolution_meters, debug=debug)
        out2 = self._detrend(out1, fit_intercept=fit_intercept, fit_dem=fit_dem, fit_coords=fit_coords,
                resolution_meters=resolution_meters, debug=debug)

        return out2

    def detrend_parallel(self, pairs=None, chunksize=None, n_jobs=-1, interactive=False, **kwargs):
        """
        Detrend and save to files a set of unwrapped interferograms combining optional topography and linear components removal
        plus optional Gaussian filtering after that.

        Parameters
        ----------
        pairs : list, tuple, array or pandas.DataFrame, optional
            A list or array of pairs of reference and repeat dates, or a DataFrame with 'ref' and 'rep' columns.
        chunksize : int or None, optional
            The number of time steps to process at a time. If None, the entire time series is processed at once.
        n_jobs : int, optional
            The number of jobs to run in parallel. -1 means using all available processors, default is -1.
        interactive : bool, optional
            Whether to return the intermediate results instead of saving them to disk. Default is False.
        **kwargs : dict
            Additional keyword arguments to be passed to the detrend function.

        Returns
        -------
        None or list
            If interactive is False (default), returns None. If interactive is True, returns a list of detrended grids.

        Examples
        --------
        Detrend plain and topography and read the results:
        sbas.detrend_parallel(pairs)
        detrended = sbas.open_grids(pairs, 'detrend')

        Detrend ionospheric effects and solid Earth's tides on large areas using Gaussian filtering
        and detrend plain and topography after that:
        sbas.detrend_parallel(pairs, wavelength=12000)

        Notes
        -----
        This function detrends and saves a set of unwrapped interferograms by combining optional topography and linear components removal.
        Additional keyword arguments can be passed to customize the detrending process.
        The detrended grids can be saved to disk or returned as intermediate results.
        """
        import dask
        import pandas as pd
        from tqdm.auto import tqdm
        import joblib

        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs = self.pairs(pairs)[['ref', 'rep']].astype(str).values

        def func(pair, **kwargs):
            #print (f'**kwargs {kwargs}')
            grid = self.open_grids([pair], 'unwrap', interactive=False)[0]
            # ignore detrending when all the fits are disabled
            out = self.detrend(grid, **kwargs)
            if interactive:
                return out
            # prepare pipeline for processing and saving
            delayed = self.save_grids([out], 'detrend', chunksize=chunksize, interactive=False)
            # perform the pipeline
            dask.persist(delayed)

        label = 'Detrending and Saving' if not interactive else 'Detrending'
        with self.tqdm_joblib(tqdm(desc=label, total=len(pairs))) as progress_bar:
            results = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(func)(pair, **kwargs) for pair in pairs)
        
        if interactive:
            return results

    def _detrend(self, dataarray, fit_intercept=True, fit_dem=True, fit_coords=True,
                resolution_meters=90, debug=False):
        """
        Detrend and return output for a single unwrapped interferogram combining optional topography and linear components removal.

        Parameters
        ----------
        dataarray : xarray.DataArray
            The input data array to detrend.
        fit_intercept : bool, optional
            Whether to remove the mean value (plane) from the data. Default is True.
        fit_dem : bool, optional
            Whether to detrend the topography. Default is True.
        fit_coords : bool, optional
            Whether to detrend the linear coordinate components. Default is True.
        resolution_meters : int, optional
            The processing resolution to prevent overfitting and reduce grid size. Default is 90.
        debug : bool, optional
            Whether to print debug information. Default is False.

        Returns
        -------
        xarray.DataArray
            The detrended 2D data array.

        Examples
        --------
        Simplest detrending:
        unwrap_detrended = sbas.detrend(pair.values[0] if isinstance(pairs, pd.DataFrame) else pair[0])

        Detrend unwrapped interferogram in radar coordinates, see for details:
        - [GitHub Issue 98](https://github.com/gmtsar/gmtsar/issues/98)
        - [GitHub Issue 411](https://github.com/gmtsar/gmtsar/issues/411)
        """
        import xarray as xr
        import numpy as np
        import dask
        from sklearn.linear_model import LinearRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        def postprocessing(out):
            return out.astype(np.float32).rename('detrend')

        # check the simplest case
        if not fit_intercept and not fit_dem and not fit_coords:
            print ('NOTE: All the detrending options disable, function does nothing')
            return dataarray

        # check simple case
        if fit_intercept and not fit_dem and not fit_coords:
            if debug:
                print ('DEBUG: Remove mean value only')
            return postprocessing(dataarray - dataarray.mean())

        # input grid can be too large
        decimator = self.pixel_decimator(resolution_meters=resolution_meters, grid=dataarray, debug=debug)
        # decimate
        dataarray_dec = decimator(dataarray)
        if debug:
            print ('DEBUG: Decimated data array', dataarray_dec.shape)

        # topography grid required to fit_dem option only
        if fit_dem:
            if debug:
                print ('DEBUG: Interpolate topography on the data grid')
            topo_ra = self.get_topo_ra()
            #topo = topo.reindex_like(unwraps[0], method='nearest')
            # use xr.zeros_like to prevent the target grid coordinates modifying
            topo = topo_ra.reindex_like(dataarray, method='nearest')
            # check chunks
            if debug:
                print ('DEBUG: regrid to resolution in meters', resolution_meters)
            # decimate
            topo_dec  = decimator(topo)
            if debug:
                print ('DEBUG: Decimated topography array', topo_dec.shape)
        else:
            topo = topo_dec = None

        # lazy calculations are useless below
        def data2fit(data, elev, yy, xx):
            y = data.values.reshape(-1) if isinstance(data, xr.DataArray) else data.reshape(-1)
            nanmask = np.isnan(y)
            # prepare regression variable
            Y = y[~nanmask]

            if fit_coords or fit_dem:
                # prepare coordinates for X regression variable
                ys = (yy.values.reshape(-1) if isinstance(yy, xr.DataArray) else yy.reshape(-1))[~nanmask]
                xs = (xx.values.reshape(-1) if isinstance(xx, xr.DataArray) else xx.reshape(-1))[~nanmask]

            if fit_dem:
                # prepare topography for X regression variable
                zs = (elev.values.reshape(-1) if isinstance(elev, xr.DataArray) else elev.reshape(-1))[~nanmask]
                zys = zs*ys
                zxs = zs*xs

            if fit_dem and fit_coords:
                X = np.column_stack([zys, zxs, ys, xs, zs])
            elif fit_dem:
                X = np.column_stack([zys, zxs, zs])
            elif fit_coords:
                X = np.column_stack([ys, xs])
            return Y, X, nanmask

        if debug:
            print ('DEBUG: linear regression calculation')

        def regr_fit():
            # build prediction model with or without plane removal (fit_intercept)
            regr = make_pipeline(StandardScaler(), LinearRegression(fit_intercept=fit_intercept))
            yy, xx = xr.broadcast(dataarray_dec.y, dataarray_dec.x)
            Y, X, _ = data2fit(dataarray_dec, topo_dec, yy, xx)

            return regr.fit(X, Y)

        # calculate for chunks
        def predict(data, elev, yy, xx, regr):
            Y, X, nanmask = data2fit(data, elev, yy, xx)
            # the chunk is NaN-filled, prediction impossible
            if nanmask.all():
                return data
            # predict when some values are not NaNs
            model = np.nan * np.zeros(data.shape)
            model.reshape(-1)[~nanmask] = regr.predict(X)
            # return data without the trend
            return data - model

        def regr_predict(regr):
            yy = xr.DataArray(dataarray.y).chunk(-1)
            xx = xr.DataArray(dataarray.x).chunk(-1)
            yy, xx = xr.broadcast(yy, xx)

            # xarray wrapper
            return xr.apply_ufunc(
                predict,
                dataarray,
                topo.chunk(dataarray.chunks) if topo is not None else None,
                yy.chunk(dataarray.chunks),
                xx.chunk(dataarray.chunks),
                dask='parallelized',
                vectorize=False,
                output_dtypes=[np.float32],
                dask_gufunc_kwargs={'regr': regr},
            )

        # build the model and return the input data without the detected trend
        return postprocessing(regr_predict(regr_fit()))

    def _gaussian(self, dataarray, wavelength, truncate=3.0, resolution_meters=90, debug=False):
        """
        Apply a lazy Gaussian filter to an input data array, allowing NaN values.

        Parameters
        ----------
        dataarray : xarray.DataArray
            The input data array with NaN values allowed.
        wavelength : float, optional
            The cut-off wavelength for the Gaussian filter in meters.
        truncate : float, optional
            The filter window size in sigmas.
        resolution_meters : float, optional
            The processing resolution for the Gaussian filter in meters.
        debug : bool, optional
            Whether to print debug information.

        Returns
        -------
        xarray.DataArray
            The filtered data array with the same coordinates as the input.

        Examples
        --------
        Detrend ionospheric effects and solid Earth's tides on a large area:
        unwrap_degaussian = sbas.degaussian(pairs.values[0], wavelength=12000)

        Lazy Gaussian filter for arrays with NaN values:
            - dataarray: Input data array with NaNs allowed.
            - wavelength: Cut-off wavelength in meters.
            - truncate: Filter window size in sigmas.
            - resolution_meters: Gaussian filter processing resolution in meters.
            - debug: Print debug information.
        Returns the filtered data array with the same coordinates as the input.
        Fast approximate calculation is silently skipped when sigma is less than 64, so the result is always exact for small filters.
        """
        import xarray as xr
        import numpy as np

        if wavelength is None:
            if debug:
                print ('DEBUG: Gaussian filter cut-off wavelength is not defined, function does nothing')
            return dataarray

        # input grid can be too large
        decimator = self.pixel_decimator(resolution_meters=resolution_meters, grid=dataarray, debug=debug)
        # decimate
        dataarray_dec = decimator(dataarray)
        # ground pixel size
        dy, dx = self.pixel_size(dataarray_dec)
        # gaussian kernel
        sigma_y = np.round(wavelength / dy)
        sigma_x = np.round(wavelength / dx)
        if debug:
            print ('DEBUG: Gaussian filtering using resolution, sigma_y, sigma_x',
                   resolution_meters, sigma_y, sigma_x)
        sigmas = (sigma_y,sigma_x)
        gaussian_dec = self.nanconvolve2d_gaussian(dataarray_dec, sigmas, truncate=truncate)
        if debug:
            print ('DEBUG: interpolate decimated filtered grid')
        gaussian = gaussian_dec.interp_like(dataarray, method='nearest')
        # revert the original chunks
        gaussian = xr.unify_chunks(dataarray, gaussian)[1]
        if debug:
            print ('DEBUG: return lazy Dask array')
        return gaussian.astype(np.float32).rename('gaussian')
