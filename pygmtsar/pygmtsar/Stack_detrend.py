# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_unwrap import Stack_unwrap

class Stack_detrend(Stack_unwrap):

    def gaussian(self, grid, wavelength, truncate=3.0, resolution_meters=60, debug=False):
        """
        Apply a lazy Gaussian filter to an input 2D or 3D data array.

        Parameters
        ----------
        dataarray : xarray.DataArray
            The input data array with NaN values allowed.
        wavelength : float
            The cut-off wavelength for the Gaussian filter in meters.
        truncate : float, optional
            Size of the Gaussian kernel, defined in terms of standard deviation, or 'sigma'. 
            It is the number of sigmas at which the window (filter) is truncated. 
            For example, if truncate = 3.0, the window will cut off at 3 sigma. Default is 3.0.
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
        Detrend ionospheric effects and solid Earth's tides on a large area and save to disk:
        stack.stack_gaussian2d(slcs, wavelength=400)
        For band-pass filtering apply the function twice and save to disk:
        model = stack.stack_gaussian2d(slcs, wavelength=400, interactive=True) \
            - stack.stack_gaussian2d(slcs, wavelength=2000, interactive=True)
        stack.save_cube(model, caption='Gaussian Band-Pass filtering')

        Detrend and return lazy xarray dataarray:
        stack.stack_gaussian2d(slcs, wavelength=400, interactive=True)
        For band-pass filtering apply the function twice:
        stack.stack_gaussian2d(slcs, wavelength=400, interactive=True) \
            - stack.stack_gaussian2d(slcs, wavelength=2000, interactive=True) 

        """
        import xarray as xr
        import numpy as np
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')

        assert self.is_ra(grid), 'ERROR: the processing requires grid in radar coordinates'
        assert np.issubdtype(grid.dtype, np.floating), 'ERROR: expected float datatype input grid'
        assert wavelength is not None, 'ERROR: Gaussian filter cut-off wavelength is not defined'

        # ground pixel size
        dy, dx = self.get_spacing(grid)
        # downscaling
        yscale, xscale = int(np.round(resolution_meters/dy)), int(np.round(resolution_meters/dx))
        # gaussian kernel
        #sigma_y = np.round(wavelength / dy / yscale, 1)
        #sigma_x = np.round(wavelength / dx / xscale, 1)
        if debug:
            print (f'DEBUG: gaussian: ground pixel size in meters: y={dy:.1f}, x={dx:.1f}')
        if (xscale <=1 and yscale <=1) or (wavelength/resolution_meters <= 3):
            # decimation is useless
            return self.multilooking(grid, wavelength=wavelength, coarsen=None, debug=debug)

        # define filter on decimated grid, the correction value is typically small
        wavelength_dec = np.sqrt(wavelength**2 - resolution_meters**2)
        if debug:
            print (f'DEBUG: gaussian: downscaling to resolution {resolution_meters}m using yscale {yscale}, xscale {xscale}')
            #print (f'DEBUG: gaussian: filtering on {resolution_meters}m grid using sigma_y0 {sigma_y}, sigma_x0 {sigma_x}')
            print (f'DEBUG: gaussian: filtering on {resolution_meters}m grid using wavelength {wavelength_dec:.1f}')

        # find stack dim
        stackvar = grid.dims[0] if len(grid.dims) == 3 else 'stack'
        #print ('stackvar', stackvar)

        # split coordinate grid to equal chunks and rest
        ys_blocks = np.array_split(grid.y, np.arange(0, grid.y.size, self.chunksize)[1:])
        xs_blocks = np.array_split(grid.x, np.arange(0, grid.x.size, self.chunksize)[1:])

        grid_dec = self.multilooking(grid, wavelength=resolution_meters, coarsen=(yscale,xscale), debug=debug)
        grid_dec_gauss = self.multilooking(grid_dec, wavelength=wavelength_dec, debug=debug)
        del grid_dec
    
        stack = []
        for stackval in grid[stackvar].values if len(grid.dims) == 3 else [None]:
            grid_in = grid_dec_gauss.sel({stackvar: stackval}) if stackval is not None else grid_dec_gauss
            grid_out = grid_in.reindex({'y': grid.y, 'x': grid.x}, method='nearest').chunk(self.chunksize)
            del grid_in
            stack.append(grid_out)
            del grid_out

        # wrap lazy Dask array to Xarray dataarray
        if len(grid.dims) == 2:
            out = stack[0]
        else:
            out = xr.concat(stack, dim=stackvar)
        del stack

        # append source grid coordinates excluding removed y, x ones
        for (k,v) in grid.coords.items():
            if k not in ['y','x']:
                out[k] = v

        return out

#     def detrend(self, dataarray, fit_intercept=True, fit_dem=True, fit_coords=True,
#                 resolution_meters=90, debug=False):
#         """
#         Detrend and return output for a single unwrapped interferogram combining optional topography and linear components removal.
# 
#         Parameters
#         ----------
#         dataarray : xarray.DataArray
#             The input data array to detrend.
#         fit_intercept : bool, optional
#             Whether to remove the mean value (plane) from the data. Default is True.
#         fit_dem : bool, optional
#             Whether to detrend the topography. Default is True.
#         fit_coords : bool, optional
#             Whether to detrend the linear coordinate components. Default is True.
#         resolution_meters : int, optional
#             The processing resolution to prevent overfitting and reduce grid size. Default is 90.
#         debug : bool, optional
#             Whether to print debug information. Default is False.
# 
#         Returns
#         -------
#         xarray.DataArray
#             The detrended 2D data array.
# 
#         Examples
#         --------
#         Simplest detrending:
#         unwrap_detrended = stack.detrend(pair.values[0] if isinstance(pairs, pd.DataFrame) else pair[0])
# 
#         Detrend unwrapped interferogram in radar coordinates, see for details:
#         - [GitHub Issue 98](https://github.com/gmtsar/gmtsar/issues/98)
#         - [GitHub Issue 411](https://github.com/gmtsar/gmtsar/issues/411)
#         """
#         import xarray as xr
#         import numpy as np
#         import dask
#         from sklearn.linear_model import LinearRegression
#         from sklearn.pipeline import make_pipeline
#         from sklearn.preprocessing import StandardScaler
# 
#         def postprocessing(out):
#             return out.astype(np.float32).rename('detrend')
# 
#         # check the simplest case
#         if not fit_intercept and not fit_dem and not fit_coords:
#             print ('NOTE: All the detrending options disable, function does nothing')
#             return dataarray
# 
#         # check simple case
#         if fit_intercept and not fit_dem and not fit_coords:
#             if debug:
#                 print ('DEBUG: Remove mean value only')
#             return postprocessing(dataarray - dataarray.mean())
# 
#         # input grid can be too large
#         decimator = self.pixel_decimator(resolution_meters=resolution_meters, grid=dataarray, debug=debug)
#         # decimate
#         dataarray_dec = decimator(dataarray)
#         if debug:
#             print ('DEBUG: Decimated data array', dataarray_dec.shape)
# 
#         # topography grid required to fit_dem option only
#         if fit_dem:
#             if debug:
#                 print ('DEBUG: Interpolate topography on the data grid')
#             topo = self.get_topo()
#             #topo = topo.reindex_like(unwraps[0], method='nearest')
#             # use xr.zeros_like to prevent the target grid coordinates modifying
#             topo = topo.reindex_like(dataarray, method='nearest')
#             # check chunks
#             if debug:
#                 print ('DEBUG: regrid to resolution in meters', resolution_meters)
#             # decimate
#             topo_dec  = decimator(topo)
#             if debug:
#                 print ('DEBUG: Decimated topography array', topo_dec.shape)
#         else:
#             topo = topo_dec = None
# 
#         # lazy calculations are useless below
#         def data2fit(data, elev, yy, xx):
#             y = data.values.reshape(-1) if isinstance(data, xr.DataArray) else data.reshape(-1)
#             nanmask = np.isnan(y)
#             # prepare regression variable
#             Y = y[~nanmask]
# 
#             if fit_coords or fit_dem:
#                 # prepare coordinates for X regression variable
#                 ys = (yy.values.reshape(-1) if isinstance(yy, xr.DataArray) else yy.reshape(-1))[~nanmask]
#                 xs = (xx.values.reshape(-1) if isinstance(xx, xr.DataArray) else xx.reshape(-1))[~nanmask]
# 
#             if fit_dem:
#                 # prepare topography for X regression variable
#                 zs = (elev.values.reshape(-1) if isinstance(elev, xr.DataArray) else elev.reshape(-1))[~nanmask]
#                 zys = zs*ys
#                 zxs = zs*xs
# 
#             if fit_dem and fit_coords:
#                 X = np.column_stack([zys, zxs, ys, xs, zs])
#             elif fit_dem:
#                 X = np.column_stack([zys, zxs, zs])
#             elif fit_coords:
#                 X = np.column_stack([ys, xs])
#             return Y, X, nanmask
# 
#         if debug:
#             print ('DEBUG: linear regression calculation')
# 
#         def regr_fit():
#             # build prediction model with or without plane removal (fit_intercept)
#             regr = make_pipeline(StandardScaler(), LinearRegression(fit_intercept=fit_intercept))
#             yy, xx = xr.broadcast(dataarray_dec.y, dataarray_dec.x)
#             Y, X, _ = data2fit(dataarray_dec, topo_dec, yy, xx)
# 
#             return regr.fit(X, Y)
# 
#         # calculate for chunks
#         def predict(data, elev, yy, xx, regr):
#             Y, X, nanmask = data2fit(data, elev, yy, xx)
#             # the chunk is NaN-filled, prediction impossible
#             if nanmask.all():
#                 return data
#             # predict when some values are not NaNs
#             model = np.nan * np.zeros(data.shape)
#             model.reshape(-1)[~nanmask] = regr.predict(X)
#             # return data without the trend
#             return data - model
# 
#         def regr_predict(regr):
#             yy = xr.DataArray(dataarray.y).chunk(-1)
#             xx = xr.DataArray(dataarray.x).chunk(-1)
#             yy, xx = xr.broadcast(yy, xx)
# 
#             # xarray wrapper
#             return xr.apply_ufunc(
#                 predict,
#                 dataarray,
#                 topo.chunk(dataarray.chunks) if topo is not None else None,
#                 yy.chunk(dataarray.chunks),
#                 xx.chunk(dataarray.chunks),
#                 dask='parallelized',
#                 vectorize=False,
#                 output_dtypes=[np.float32],
#                 dask_gufunc_kwargs={'regr': regr},
#             )
# 
#         # build the model and return the input data without the detected trend
#         return postprocessing(regr_predict(regr_fit()))

#     def gaussian(self, grid, wavelength, truncate=3.0, resolution_meters=90, debug=False):
#         """
#         Apply a lazy Gaussian filter to an input 2D or 3D data array.
# 
#         Parameters
#         ----------
#         dataarray : xarray.DataArray
#             The input data array with NaN values allowed.
#         wavelength : float
#             The cut-off wavelength for the Gaussian filter in meters.
#         truncate : float, optional
#             Size of the Gaussian kernel, defined in terms of standard deviation, or 'sigma'. 
#             It is the number of sigmas at which the window (filter) is truncated. 
#             For example, if truncate = 3.0, the window will cut off at 3 sigma. Default is 3.0.
#         resolution_meters : float, optional
#             The processing resolution for the Gaussian filter in meters.
#         debug : bool, optional
#             Whether to print debug information.
# 
#         Returns
#         -------
#         xarray.DataArray
#             The filtered data array with the same coordinates as the input.
# 
#         Examples
#         --------
#         Detrend ionospheric effects and solid Earth's tides on a large area and save to disk:
#         stack.stack_gaussian2d(slcs, wavelength=400)
#         For band-pass filtering apply the function twice and save to disk:
#         model = stack.stack_gaussian2d(slcs, wavelength=400, interactive=True) \
#             - stack.stack_gaussian2d(slcs, wavelength=2000, interactive=True)
#         stack.save_cube(model, caption='Gaussian Band-Pass filtering')
# 
#         Detrend and return lazy xarray dataarray:
#         stack.stack_gaussian2d(slcs, wavelength=400, interactive=True)
#         For band-pass filtering apply the function twice:
#         stack.stack_gaussian2d(slcs, wavelength=400, interactive=True) \
#             - stack.stack_gaussian2d(slcs, wavelength=2000, interactive=True) 
# 
#         """
#         import xarray as xr
#         import numpy as np
# 
#         assert self.is_ra(grid), 'ERROR: the processing requires grid in radar coordinates'
#         assert wavelength is not None, 'ERROR: Gaussian filter cut-off wavelength is not defined'
# 
#         # ground pixel size
#         dy, dx = self.get_spacing(grid)
#         # reduction
#         yscale, xscale = int(np.round(resolution_meters/dy)), int(np.round(resolution_meters/dx))
#         # gaussian kernel
#         sigma_y = np.round(wavelength / dy / yscale)
#         sigma_x = np.round(wavelength / dx / xscale)
#         if debug:
#             print (f'DEBUG: average ground pixel size in meters: y={dy}, x={dx}')
#             print (f'DEBUG: yscale {yscale}, xscale {xscale} to resolution {resolution_meters} m')
#             print ('DEBUG: Gaussian filtering using resolution, sigma_y, sigma_x', resolution_meters, sigma_y, sigma_x)
# 
#         # find stack dim
#         stackvar = grid.dims[0] if len(grid.dims) == 3 else 'stack'
#         #print ('stackvar', stackvar)
#     
#         stack = []
#         for stackval in grid[stackvar].values if len(grid.dims) == 3 else [None]:
#             block = grid.sel({stackvar: stackval}) if stackval is not None else grid
#             block_dec = self.antialiasing_downscale(block, wavelength=resolution_meters, coarsen=(yscale,xscale), debug=debug)
#             gaussian_dec = self.nanconvolve2d_gaussian(block_dec, (sigma_y,sigma_x), truncate=truncate)
#             # interpolate decimated filtered grid to original resolution
#             gaussian = gaussian_dec.interp_like(block, method='nearest')
#             # revert the original chunks
#             gaussian = xr.unify_chunks(block, gaussian)[1]
#             stack.append(gaussian.astype(np.float32).rename('gaussian'))
# 
#         # wrap lazy Dask array to Xarray dataarray
#         if len(grid.dims) == 2:
#             out = stack[0]
#         else:
#             out = xr.concat(stack, dim=stackvar)
#         del stack
# 
#         # append source grid coordinates excluding removed y, x ones
#         for (k,v) in grid.coords.items():
#             if k not in ['y','x']:
#                 out[k] = v
#     
#         return out
