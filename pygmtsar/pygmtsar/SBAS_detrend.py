#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_unwrap import SBAS_unwrap

class SBAS_detrend(SBAS_unwrap):

    def detrend_parallel(self, pairs, n_jobs=-1, **kwargs):
        from tqdm.auto import tqdm
        import joblib

        def func(pair, **kwargs):
            #print (f'**kwargs {kwargs}')
            grid = self.open_grids([pair], 'unwrap', interactive=False)
            
            kwargs1 = {k:v for k,v in kwargs.items() if k in ['wavelength', 'truncate', 'resolution_meters', 'debug']}
            #print (f'kwargs1 {kwargs1}')
            if 'wavelength' in kwargs1:
                out1 = self.degaussian(grid[0], **kwargs1)
            else:
                out1 = grid[0]
        
            kwargs2 = {k:v for k,v in kwargs.items() if k in ['fit_intercept', 'fit_dem', 'fit_coords',
                                                              'resolution_meters', 'debug']}
            #print (f'kwargs2 {kwargs2}')
            out2 = self.detrend(out1, **kwargs2)
        
            self.save_grids([out2], 'detrend', interactive=False)

        with self.tqdm_joblib(tqdm(desc='Detrending and Saving', total=len(pairs))) as progress_bar:
            joblib.Parallel(n_jobs=1)(joblib.delayed(func)(pair, **kwargs) for pair in pairs.values)

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
        assert fit_intercept or fit_dem or fit_coords, 'All the detrending options disable, function does nothing'

        # check simple case
        if fit_intercept and not fit_dem and not fit_coords:
            if debug:
                print ('DEBUG: Remove mean value only')
            return postprocessing(dataarray - dataarray.mean())

        # input grid can be too large
        decimator = self.pixel_decimator(resolution_meters=resolution_meters, grid=dataarray, debug=debug)
        # decimate
        dataarray_dec = decimator(dataarray)

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
        else:
            topo = topo_dec = None

        # lazy calculations are useless below
        def data2fit(data, grid):
            y = data.values.reshape(-1)
            nanmask = np.isnan(y)
            # prepare regression variable
            Y = y[~nanmask]

            if fit_coords or fit_dem:
                # prepare coordinates for X regression variable
                yy, xx = xr.broadcast(data.y, data.x)
                ys = yy.values.reshape(-1)[~nanmask]
                xs = xx.values.reshape(-1)[~nanmask]

            if fit_dem:
                # prepare topography for X regression variable
                zs = grid.values.reshape(-1)[~nanmask]
                zys = zs*ys
                zxs = zs*xs

            if fit_dem and fit_coords:
                if debug:
                    print ('DEBUG: Detrend topography and coordinates')
                X = np.column_stack([zys, zxs, ys, xs, zs])
            elif fit_dem:
                if debug:
                    print ('DEBUG: Detrend topography only')
                X = np.column_stack([zys, zxs, zs])
            elif fit_coords:
                if debug:
                    print ('DEBUG: Detrend coordinates only')
                X = np.column_stack([ys, xs])
            return Y, X, nanmask

        # build prediction model with or without plane removal (fit_intercept)
        regr = make_pipeline(StandardScaler(), LinearRegression(fit_intercept=fit_intercept))
        Y, X, _ = data2fit(dataarray_dec, topo_dec)
        regr.fit(X, Y)
    
        # TODO: calculate for chunks
        Y, X, nanmask = data2fit(dataarray, topo)
        model = xr.full_like(dataarray, np.nan).compute()
        if debug:
            print ('DEBUG: model', model)
        model.data.reshape(-1)[~nanmask] = regr.predict(X)
        return postprocessing(dataarray - model)

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
