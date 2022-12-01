#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_unwrap import SBAS_unwrap

class SBAS_detrend(SBAS_unwrap):

    def detrend_parallel(self, pairs=None, n_jobs=-1, interactive=False, **kwargs):
        import dask
        import pandas as pd
        from tqdm.auto import tqdm
        import joblib

        if pairs is None:
            pairs = self.find_pairs()
        elif isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        def func(pair, **kwargs):
            #print (f'**kwargs {kwargs}')
            grid = self.open_grids([pair], 'unwrap', interactive=False)

            # without any processing options return the input as is
            out = grid[0]

            if 'wavelength' in kwargs:
                kwargs1 = {k:v for k,v in kwargs.items() if k in ['wavelength', 'truncate', 'resolution_meters', 'debug']}
                #print (f'kwargs1 {kwargs1}')
                out = self.degaussian(out, **kwargs1)

            kwargs2 = {k:v for k,v in kwargs.items() if k in ['fit_intercept', 'fit_dem', 'fit_coords',
                                                              'resolution_meters', 'debug']}
            #print (f'kwargs2 {kwargs2}')
            # ignore detrending when all the fits are disabled
            out = self.detrend(out, **kwargs2)

            if interactive:
                return out

            # prepare pipeline for processing and saving
            delayed = self.save_grids([out], 'detrend', interactive=False)

            # perform the pipeline
            dask.persist(delayed)

        label = 'Detrending and Saving' if not interactive else 'Detrending'
        with self.tqdm_joblib(tqdm(desc=label, total=len(pairs))) as progress_bar:
            results = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(func)(pair, **kwargs) for pair in pairs)
        
        if interactive:
            return results

    def detrend(self, dataarray, fit_intercept=True, fit_dem=True, fit_coords=True,
                resolution_meters=90, debug=False):
        """
        Detrend unwrapped interferogram in radar coordinates, see for details
            https://github.com/gmtsar/gmtsar/issues/98
            https://github.com/gmtsar/gmtsar/issues/411
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
            y = data.reshape(-1)
            nanmask = np.isnan(y)
            # prepare regression variable
            Y = y[~nanmask]

            if fit_coords or fit_dem:
                # prepare coordinates for X regression variable
                ys = yy.reshape(-1)[~nanmask]
                xs = xx.reshape(-1)[~nanmask]

            if fit_dem:
                # prepare topography for X regression variable
                zs = elev.reshape(-1)[~nanmask]
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
            Y, X, _ = data2fit(dataarray_dec.values, topo_dec.values, yy.values, xx.values)
        
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
                topo.chunk(dataarray.chunks),
                yy.chunk(dataarray.chunks),
                xx.chunk(dataarray.chunks),
                dask='parallelized',
                vectorize=False,
                output_dtypes=[np.float32],
                dask_gufunc_kwargs={'regr': regr},
            )
    
        # build the model and return the input data without the detected trend
        return postprocessing(regr_predict(regr_fit()))

    def degaussian(self, dataarray, wavelength, truncate=3.0, resolution_meters=90, debug=False):
        """
        Lazy Gaussian filter for arrays with NaN values.
            dataarray - input dataarray with NaNs allowed,
            wavelength - cut-off wavelength [m],
            truncate - filter window size [sigma],
            resolution_meters - Gaussian filter processing resolution [m],
            debug - print debug information.
        Returns filtered dataarray with the same coordinates as input one.
        Fast approximate calculation silently skipped when sigma is less than 64 so the result is always exact for small filters.
        """
        import xarray as xr
        import numpy as np

        if wavelength is None:
            print ('NOTE: Gaussian filter cut-off wavelength is not defined, function does nothing')
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
        return (dataarray - gaussian).astype(np.float32).rename('degaussian')
