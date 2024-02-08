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
# 
#     @staticmethod
#     def regression_linear(data, grid, fit_intercept=True):
#         import numpy as np
#         import xarray as xr
#         import dask
# 
#         # define topography on the same grid and fill NaNs
#         grid = grid.reindex_like(data, method='nearest').fillna(0)
# 
#         # find stack dim
#         stackvar = data.dims[0] if len(data.dims) == 3 else 'stack'
#         #print ('stackvar', stackvar)
#         shape2d = data.shape[1:] if len(data.dims) == 3 else data.shape
#         #print ('shape2d', shape2d)
# 
#         @dask.delayed
#         def regression_block(stackval, fit_intercept):
#             from sklearn.linear_model import LinearRegression
#             from sklearn.pipeline import make_pipeline
#             from sklearn.preprocessing import StandardScaler
# 
#             # use outer variable
#             data_block = (data.sel({stackvar: stackval}) if stackval is not None else data).compute(n_workers=1)
# 
#             grid_values = grid.values.ravel()
#             data_values = data_block.values.ravel()
#             nanmask = np.isnan(grid_values) | np.isnan(data_values)
#         
#             # build prediction model with or without plane removal (fit_intercept)
#             regr = make_pipeline(StandardScaler(), LinearRegression(fit_intercept=fit_intercept))
#             # fit 1D non-NaNs phase and topography and predict on the non-NaNs 2D topography grid
#             data_grid = regr.fit(np.column_stack([grid_values [~nanmask]]),
#                                  np.column_stack([data_values[~nanmask]]))\
#                         .predict(np.column_stack([grid_values])).reshape(data_block.shape)
#             # cleanup
#             del grid_values, data_values, nanmask, data_block
#             return data_grid
#         
#         stack = []
#         for stackval in data[stackvar].values if len(data.dims) == 3 else [None]:
#             #print ('stackval', stackval)
#             block = dask.array.from_delayed(regression_block(stackval, fit_intercept=fit_intercept),
#                                             shape=shape2d, dtype=np.float32)
#             stack.append(block)
#             del block
# 
#         return xr.DataArray(dask.array.stack(stack) if len(data.dims) == 3 else stack[0],
#                             coords=data.coords)\
#                .rename(data.name)

#     @staticmethod
#     def regression_linear(data, variables, weight=None, valid_pixels_threshold=10000, fit_intercept=True):
#         """
#         topo = sbas.get_topo().coarsen({'x': 4}, boundary='trim').mean()
#         yy, xx = xr.broadcast(topo.y, topo.x)
#         strat_sbas = sbas.regression_linear(unwrap_sbas.phase,
#                 [topo,    topo*yy,    topo*xx,    topo*yy*xx,
#                  topo**2, topo**2*yy, topo**2*xx, topo**2*yy*xx,
#                  topo**3, topo**3*yy, topo**3*xx, topo**3*yy*xx,
#                  yy, xx,
#                  yy**2, xx**2, yy*xx,
#                  yy**3, xx**3, yy**2*xx, xx**2*yy], corr_sbas)
#         """
#         import numpy as np
#         import xarray as xr
#         import dask
#         from sklearn.linear_model import LinearRegression
#         from sklearn.pipeline import make_pipeline
#         from sklearn.preprocessing import StandardScaler
# 
#         if weight is not None and len(weight.dims) == 3 and weight.shape != data.shape:
#             raise Exception(f'Argument "weight" 3D shape {weight.shape} should be equal to "data" 3D shape {data.shape}')
#         if weight is not None and len(weight.dims) == 2 and weight.shape != data.shape[1:]:
#             raise Exception(f'Argument "weight" 2D shape {weight.shape} should be equal to "data" 2D shape {data.shape[1:]}')
# 
#         # find stack dim
#         stackvar = data.dims[0] if len(data.dims) == 3 else 'stack'
#         #print ('stackvar', stackvar)
#         shape2d = data.shape[1:] if len(data.dims) == 3 else data.shape
#         #print ('shape2d', shape2d)
#         chunk2d = data.chunks[1:] if len(data.dims) == 3 else data.chunks
#         #print ('chunk2d', chunk2d)
#     
#         if isinstance(variables, (list, tuple)):
#             variables = xr.concat(variables, dim='stack')
#         elif not isinstance(variables, xr.DataArray) or len(variables.dims) not in (2, 3):
#             raise Exception('Argument "variables" should be 2D or 3D Xarray dataarray of list of 2D Xarray dataarrays')
#         elif len(variables.dims) == 2:
#             variables = variables.expand_dims('stack')
#         elif len(variables.dims) == 3 and not variables.dims[0] == 'stack':
#             raise Exception('Argument "variables" 3D Xarray dataarray needs the first dimension name "stack"')
#         #print ('variables', variables)
# 
#         def regression_block(data, variables, weight, valid_pixels_threshold, fit_intercept):
#             data_values  = data.ravel()
#             variables_values = variables.reshape(-1, variables.shape[-1]).T
#             #assert 0, f'TEST: {data_values.shape}, {variables_values.shape}, {weight.shape}'
#             if weight.size > 1:
#                 weight_values = weight.ravel()
#                 nanmask = np.isnan(data_values) | np.isnan(weight_values) | np.any(np.isnan(variables_values), axis=0)
#             else:
#                 weight_values = None
#                 nanmask = np.isnan(data_values) | np.any(np.isnan(variables_values), axis=0)
# 
#             # regression requires enough amount of valid pixels
#             if data_values.size - np.sum(nanmask) < valid_pixels_threshold:
#                 del data_values, variables_values, weight_values, nanmask
#                 return np.nan * np.zeros(data.shape)
# 
#             # build prediction model with or without plane removal (fit_intercept)
#             regr = make_pipeline(StandardScaler(), LinearRegression(fit_intercept=fit_intercept, copy_X=False, n_jobs=1))
#             fit_params = {'linearregression__sample_weight': weight_values[~nanmask]} if weight.size > 1 else {}
#             regr.fit(variables_values[:, ~nanmask].T, data_values[~nanmask], **fit_params)
#             del weight_values, data_values
#             model = np.full_like(data, np.nan).ravel()
#             model[~nanmask] = regr.predict(variables_values[:, ~nanmask].T)
#             del variables_values, regr
#             return model.reshape(data.shape)
# 
#         # xarray wrapper
#         model = xr.apply_ufunc(
#             regression_block,
#             data,
#             variables.chunk(dict(stack=-1, y=chunk2d[0], x=chunk2d[1])),
#             weight.chunk(dict(y=chunk2d[0], x=chunk2d[1])) if weight is not None else weight,
#             dask='parallelized',
#             vectorize=False,
#             output_dtypes=[np.float32],
#             input_core_dims=[[], ['stack'], []],
#             dask_gufunc_kwargs={'valid_pixels_threshold': valid_pixels_threshold, 'fit_intercept': fit_intercept},
#         )
#         return model

    @staticmethod
    def regression(data, variables, weight=None, valid_pixels_threshold=1000, **kwargs):
        """
        topo = sbas.get_topo().coarsen({'x': 4}, boundary='trim').mean()
        yy, xx = xr.broadcast(topo.y, topo.x)
        strat_sbas = sbas.regression(unwrap_sbas.phase,
                [topo,    topo*yy,    topo*xx,    topo*yy*xx,
                 topo**2, topo**2*yy, topo**2*xx, topo**2*yy*xx,
                 topo**3, topo**3*yy, topo**3*xx, topo**3*yy*xx,
                 yy, xx,
                 yy**2, xx**2, yy*xx,
                 yy**3, xx**3, yy**2*xx, xx**2*yy], corr_sbas)


        topo = sbas.interferogram(topophase)
        inc = decimator(sbas.incidence_angle())
        yy, xx = xr.broadcast(topo.y, topo.x)
        variables = [topo,  topo*yy,  topo*xx, topo*yy*xx]
        trend_sbas = sbas.regression(intf_sbas, variables, corr_sbas)
        """
        import numpy as np
        import xarray as xr
        import dask
        # 'linear'
        from sklearn.linear_model import LinearRegression
        # 'sgd'
        from sklearn.linear_model import SGDRegressor
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        # find stack dim
        stackvar = data.dims[0] if len(data.dims) >= 3 else 'stack'
        #print ('stackvar', stackvar)
        shape2d = data.shape[1:] if len(data.dims) == 3 else data.shape
        #print ('shape2d', shape2d)
        chunk2d = data.chunks[1:] if len(data.dims) == 3 else data.chunks
        #print ('chunk2d', chunk2d)

        def regression_block(data, weight, *args, **kwargs):
            #assert 0, stack.shape
            data_values  = data.ravel().astype(np.float64)
            # manage variable number of variables
            variables_stack = np.stack([arg[0] if arg.ndim==3 else arg for arg in args])
            #variables_values = variables_stack.reshape(-1, variables_stack.shape[0]).T
            variables_values = variables_stack.reshape(variables_stack.shape[0], -1)
            del variables_stack
            #assert 0, f'TEST: {data_values.shape}, {variables_values.shape}'

            nanmask_data = ~np.isfinite(data_values)
            nanmask_values = np.any(~np.isfinite(variables_values), axis=0)
            if weight.size > 1:
                weight_values = weight.ravel().astype(np.float64)
                nanmask_weight = ~np.isfinite(weight_values)
                nanmask = nanmask_data | nanmask_values | nanmask_weight
                #assert 0, f'TEST weight: {data_values.shape}, {variables_values.shape}, {weight_values.shape}'
            else:
                weight_values = None
                nanmask_weight = None
                nanmask = nanmask_data | nanmask_values

            # regression requires enough amount of valid pixels
            if data_values.size - np.sum(nanmask) < valid_pixels_threshold:
                del data_values, variables_values, weight_values
                del nanmask_data, nanmask_values, nanmask_weight, nanmask
                return np.nan * np.zeros(data.shape)

            # build prediction model with or without plane removal (fit_intercept)
            algorithm = kwargs.pop('algorithm', 'linear')
            if algorithm == 'sgd':
                regr = make_pipeline(StandardScaler(), SGDRegressor(**kwargs))
                fit_params = {'sgdregressor__sample_weight': weight_values[~nanmask]} if weight.size > 1 else {}
            elif algorithm == 'linear':
                regr = make_pipeline(StandardScaler(), LinearRegression(**kwargs, copy_X=False, n_jobs=1))
                fit_params = {'linearregression__sample_weight': weight_values[~nanmask]} if weight.size > 1 else {}
            else:
                raise ValueError(f"Unsupported algorithm {kwargs['algorithm']}. Should be 'linear' or 'sgd'")
            regr.fit(variables_values[:, ~nanmask].T, data_values[~nanmask], **fit_params)
            del weight_values, data_values
            model = np.full_like(data, np.nan).ravel()
            model[~nanmask_values] = regr.predict(variables_values[:, ~nanmask_values].T)
            del variables_values, regr
            del nanmask_data, nanmask_values, nanmask_weight, nanmask
            return model.reshape(data.shape).astype(np.float32)

        dshape = data[0].shape if data.ndim==3 else data.shape
        if isinstance(variables, (list, tuple)):
            vshapes = [v[0].shape if v.ndim==3 else v.shape for v in variables]
            equals = np.all([vshape == dshape for vshape in vshapes])
            assert equals, f'ERROR: shapes of variables slices {vshapes} and data slice {dshape} differ.'
            #assert equals, f'{dshape} {vshapes}, {equals}'
            variables_stack = [v.chunk(dict(y=chunk2d[0], x=chunk2d[1])) for v in variables]
        else:
            vshape = variables[0].shape if variables.ndim==3 else variables.shape
            assert {vshape} == {dshape}, f'ERROR: shapes of variables slice {vshape} and data slice {dshape} differ.'
            variables_stack = [variables.chunk(dict(y=chunk2d[0], x=chunk2d[1]))]

        # xarray wrapper
        model = xr.apply_ufunc(
            regression_block,
            data,
            weight.chunk(dict(y=chunk2d[0], x=chunk2d[1])) if weight is not None else weight,
            *variables_stack,
            dask='parallelized',
            vectorize=False,
            output_dtypes=[np.float32],
            dask_gufunc_kwargs={**kwargs},
        )
        del variables_stack
        return model

    def regression_linear(self, data, variables, weight=None, valid_pixels_threshold=1000, fit_intercept=True):
        """   
        topo = sbas.get_topo().coarsen({'x': 4}, boundary='trim').mean()
        yy, xx = xr.broadcast(topo.y, topo.x)
        strat_sbas = sbas.regression_linear(unwrap_sbas.phase,
                [topo,    topo*yy,    topo*xx,    topo*yy*xx,
                 topo**2, topo**2*yy, topo**2*xx, topo**2*yy*xx,
                 topo**3, topo**3*yy, topo**3*xx, topo**3*yy*xx,
                 yy, xx,
                 yy**2, xx**2, yy*xx,
                 yy**3, xx**3, yy**2*xx, xx**2*yy], corr_sbas)
        """
        return self.regression(data, variables, weight, valid_pixels_threshold, algorithm='linear',
                                fit_intercept=fit_intercept)

    def regression_sgd(self, data, variables, weight=None, valid_pixels_threshold=1000,
                      max_iter=1000, tol=1e-3, alpha=0.0001, l1_ratio=0.15):
        """
        Perform regression on a dataset using the SGDRegressor model from scikit-learn.

        This function applies Stochastic Gradient Descent (SGD) regression to fit the given 'data' against a set of 'variables'. It's suitable for large datasets and handles high-dimensional features efficiently.

        Parameters:
        data (xarray.DataArray): The target data array to fit.
        variables (xarray.DataArray or list of xarray.DataArray): Predictor variables. It can be a single 3D DataArray or a list of 2D DataArrays.
        weight (xarray.DataArray, optional): Weights for each data point. Useful if certain data points are more important. Defaults to None.
        valid_pixels_threshold (int, optional): Minimum number of valid pixels required for the regression to be performed. Defaults to 10000.
        max_iter (int, optional): Maximum number of passes over the training data (epochs). Defaults to 1000.
        tol (float, optional): Stopping criterion. If not None, iterations will stop when (loss > best_loss - tol) for n_iter_no_change consecutive epochs. Defaults to 1e-3.
        alpha (float, optional): Constant that multiplies the regularization term. Higher values mean stronger regularization. Defaults to 0.0001.
        l1_ratio (float, optional): The Elastic Net mixing parameter, with 0 <= l1_ratio <= 1. l1_ratio=0 corresponds to L2 penalty, l1_ratio=1 to L1. Defaults to 0.15.

        Returns:
        xarray.DataArray: The predicted values as an xarray DataArray, fitted by the SGDRegressor model.

        Notes:
        - SGDRegressor is well-suited for large datasets due to its efficiency in handling large-scale and high-dimensional data.
        - Proper tuning of parameters (max_iter, tol, alpha, l1_ratio) is crucial for optimal performance and prevention of overfitting.

        Example:  
        decimator = sbas.decimator(resolution=15, grid=(1,1))
        topo = decimator(sbas.get_topo())
        inc = decimator(sbas.incidence_angle())
        yy, xx = xr.broadcast(topo.y, topo.x)
        trend_sbas = sbas.regression(unwrap_sbas.phase,
                [topo,    topo*yy,    topo*xx,    topo*yy*xx,
                 topo**2, topo**2*yy, topo**2*xx, topo**2*yy*xx,
                 topo**3, topo**3*yy, topo**3*xx, topo**3*yy*xx,
                 inc,     inc**yy,    inc*xx,     inc*yy*xx,
                 yy, xx,
                 yy**2, xx**2, yy*xx,
                 yy**3, xx**3, yy**2*xx, xx**2*yy], corr_sbas)
        """
        return self.regression(data, variables, weight, valid_pixels_threshold, algorithm='sgd',
                               max_iter=max_iter, tol=tol, alpha=alpha, l1_ratio=l1_ratio)

    def polyfit(self, data, weight=None, degree=0, days=None, count=None):
        import xarray as xr
        import pandas as pd
        import numpy as np

        pairs, dates = self.get_pairs(data, dates=True)

        models = []
        for date in dates:
            #print ('date', date)
            data_pairs = pairs[(pairs.ref==date)|(pairs.rep==date)].pair.values
            #print ('data_pairs', data_pairs)
            if weight is None:
                stack = data.sel(pair=data_pairs)
            else:
                stack = data.sel(pair=data_pairs) * np.sqrt(weight.sel(pair=data_pairs))
            stack_days = xr.where(stack.ref < pd.Timestamp(date),
                           (stack.ref - stack.rep).dt.days,
                           (stack.rep - stack.ref).dt.days)
            #print ('stack_days', stack_days)
            # select smallest intervals
            stack_days_selected = stack_days[np.argsort(np.abs(stack_days.values))][:count]
            if days is not None:
                stack_days_selected = stack_days_selected[np.abs(stack_days_selected)<=days]
            #print ('stack_days_selected', stack_days_selected)
            linear_fit = (np.sign(stack_days)*stack).assign_coords(time=stack_days)\
                [stack.pair.isin(stack_days_selected.pair)]\
                .swap_dims({'pair': 'time'})\
                .sortby(['ref', 'rep'])\
                .polyfit(dim='time', deg=degree)
            model = linear_fit.polyfit_coefficients.sel(degree=degree)
            models.append(model.assign_coords(date=pd.to_datetime(date)))
            del data_pairs, stack, stack_days, stack_days_selected, linear_fit, model
        model = xr.concat(models, dim='date')
        del models

        out = xr.concat([(model.sel(date=ref).drop('date') - model.sel(date=rep).drop('date'))\
                                 .assign_coords(pair=str(ref.date()) + ' ' + str(rep.date()), ref=ref, rep=rep) \
                          for ref, rep in zip(pairs['ref'], pairs['rep'])], dim='pair')
        return out.rename(data.name)

#     def polyfit(self, data, weight=None, degree=0, variable=None, count=None):
#         import xarray as xr
#         import pandas as pd
#         import numpy as np
# 
#         pairs, dates = self.get_pairs(data, dates=True)
#         if variable is not None:
#             pairs['variable'] = variable.values if isinstance(variable, pd.Series) else variable
# 
#         models = []
#         for date in dates:
#             data_pairs = pairs[(pairs.ref==date)|(pairs.rep==date)].pair.values
#             #print (data_pairs)
#             stack = data.sel(pair=data_pairs)
#             if weight is None:
#                 stack = data.sel(pair=data_pairs)
#             else:
#                 stack = data.sel(pair=data_pairs) * np.sqrt(weight.sel(pair=data_pairs))
#             # days, positive and negative
#             days = xr.where(stack.ref < pd.Timestamp(date),
#                            (stack.ref - stack.rep).dt.days,
#                            (stack.rep - stack.ref).dt.days)
#             # select smallest intervals
#             days_selected = days[np.argsort(np.abs(days.values))][:count]
#             if variable is not None:
#                 data_variables = xr.DataArray(pairs[pairs.pair.isin(data_pairs)]['variable'].values, coords=days.coords)
#             else:
#                 data_variables = days
#             #print (days, days_selected, data_variables)
#             #.sortby(['ref', 'rep'])
#             linear_fit = (np.sign(days)*stack)\
#                 .assign_coords(time=data_variables)[stack.pair.isin(days_selected.pair)]\
#                 .swap_dims({'pair': 'variable'})\
#                 .polyfit(dim='variable', deg=degree)
#             model = linear_fit.polyfit_coefficients.sel(degree=degree)
#             models.append(model.assign_coords(date=pd.to_datetime(date)))
#             del data_pairs, stack, days, days_selected, data_variables, linear_fit, model
#         model = xr.concat(models, dim='date')
#         del models
# 
#         out = xr.concat([(model.sel(date=ref).drop('date') - model.sel(date=rep).drop('date'))\
#                                  .assign_coords(pair=str(ref.date()) + ' ' + str(rep.date()), ref=ref, rep=rep) \
#                           for ref, rep in zip(pairs['ref'], pairs['rep'])], dim='pair')
# 
#         return out.rename(data.name)

    def gaussian(self, grid, wavelength, truncate=3.0, resolution=60, debug=False):
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
        resolution : float, optional
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
#         import warnings
#         # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
#         warnings.filterwarnings('ignore')
#         warnings.filterwarnings('ignore', module='dask')
#         warnings.filterwarnings('ignore', module='dask.core')

        assert self.is_ra(grid), 'ERROR: the processing requires grid in radar coordinates'
        assert np.issubdtype(grid.dtype, np.floating), 'ERROR: expected float datatype input grid'
        assert wavelength is not None, 'ERROR: Gaussian filter cut-off wavelength is not defined'

        # ground pixel size
        dy, dx = self.get_spacing(grid)
        # downscaling
        yscale, xscale = int(np.round(resolution/dy)), int(np.round(resolution/dx))
        # gaussian kernel
        #sigma_y = np.round(wavelength / dy / yscale, 1)
        #sigma_x = np.round(wavelength / dx / xscale, 1)
        if debug:
            print (f'DEBUG: gaussian: ground pixel size in meters: y={dy:.1f}, x={dx:.1f}')
        if (xscale <=1 and yscale <=1) or (wavelength/resolution <= 3):
            # decimation is useless
            return self.multilooking(grid, wavelength=wavelength, coarsen=None, debug=debug)

        # define filter on decimated grid, the correction value is typically small
        wavelength_dec = np.sqrt(wavelength**2 - resolution**2)
        if debug:
            print (f'DEBUG: gaussian: downscaling to resolution {resolution}m using yscale {yscale}, xscale {xscale}')
            #print (f'DEBUG: gaussian: filtering on {resolution}m grid using sigma_y0 {sigma_y}, sigma_x0 {sigma_x}')
            print (f'DEBUG: gaussian: filtering on {resolution}m grid using wavelength {wavelength_dec:.1f}')

        # find stack dim
        stackvar = grid.dims[0] if len(grid.dims) == 3 else 'stack'
        #print ('stackvar', stackvar)

        # split coordinate grid to equal chunks and rest
        ys_blocks = np.array_split(grid.y, np.arange(0, grid.y.size, self.chunksize)[1:])
        xs_blocks = np.array_split(grid.x, np.arange(0, grid.x.size, self.chunksize)[1:])

        grid_dec = self.multilooking(grid, wavelength=resolution, coarsen=(yscale,xscale), debug=debug)
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
#                 resolution=90, debug=False):
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
#         resolution : int, optional
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
#         decimator = self.pixel_decimator(resolution=resolution, grid=dataarray, debug=debug)
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
#                 print ('DEBUG: regrid to resolution in meters', resolution)
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

#     def gaussian(self, grid, wavelength, truncate=3.0, resolution=90, debug=False):
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
#         resolution : float, optional
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
#         yscale, xscale = int(np.round(resolution/dy)), int(np.round(resolution/dx))
#         # gaussian kernel
#         sigma_y = np.round(wavelength / dy / yscale)
#         sigma_x = np.round(wavelength / dx / xscale)
#         if debug:
#             print (f'DEBUG: average ground pixel size in meters: y={dy}, x={dx}')
#             print (f'DEBUG: yscale {yscale}, xscale {xscale} to resolution {resolution} m')
#             print ('DEBUG: Gaussian filtering using resolution, sigma_y, sigma_x', resolution, sigma_y, sigma_x)
# 
#         # find stack dim
#         stackvar = grid.dims[0] if len(grid.dims) == 3 else 'stack'
#         #print ('stackvar', stackvar)
#     
#         stack = []
#         for stackval in grid[stackvar].values if len(grid.dims) == 3 else [None]:
#             block = grid.sel({stackvar: stackval}) if stackval is not None else grid
#             block_dec = self.antialiasing_downscale(block, wavelength=resolution, coarsen=(yscale,xscale), debug=debug)
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

#     def turbulence(self, phase, date_crop=0, symmetrical=False):
#         pairs, dates = self.get_pairs(phase, dates=True)
#     
#         turbos = []
#         for date in dates:
#             ref = pairs[pairs.ref==date]
#             rep = pairs[pairs.rep==date]
#             # calculate left and right pairs
#             if symmetrical:
#                 count = min(len(ref), len(rep))
#             else:
#                 count = None
#             #print (date, len(ref), len(rep), '=>', count)
#             if len(ref) < 1 or len(rep) < 1:
#                 # correction calculation is not possible
#                 #turbo = xr.zeros_like((detrend - phase_topo).isel(pair=0)).drop(['pair','ref','rep'])
#                 continue
#             else:
#                 ref_data = phase.sel(pair=ref.pair.values).isel(pair=slice(None,count))
#                 #ref_weight = 1/corr60m.sel(pair=ref_data.pair)
#                 #print (ref_data)
#                 rep_data = phase.sel(pair=rep.pair.values ).isel(pair=slice(None,count))
#                 #rep_weight = 1/corr60m.sel(pair=rep_data.pair)
#                 #print (rep_data)
#                 #mask = (ref_weight.sum('pair') + rep_weight.sum('pair'))/2/count
#                 turbo = (ref_data.mean('pair') - rep_data.mean('pair')) / 2
#                 #turbo = ((ref_data*ref_weight).mean('pair')/ref_weight.sum('pair') -\
#                 #         (rep_data*rep_weight).mean('pair')/rep_weight.sum('pair')) / 2
#                 #.where(mask>MASK)
#             turbos.append(turbo.assign_coords(date=pd.to_datetime(date)))
#             del turbo
#         turbo = xr.concat(turbos, dim='date')
#         del turbos
#     
#         # convert dates to pairs
#         #pairs = pairs[pairs.ref.isin(turbo.date.values)&pairs.rep.isin(turbo.date.values)]
#         fake = xr.zeros_like(phase.isel(pair=0))
# 
#     #     return [(
#     #         (turbo.sel(date=ref).drop('date') if ref in turbo.date else fake) - \
#     #         (turbo.sel(date=rep).drop('date') if rep in turbo.date else fake)
#     #     ).assign_coords(pair=str(ref.date()) + ' ' + str(rep.date()), ref=ref, rep=rep) \
#     #     for ref, rep in zip(pairs['ref'], pairs['rep'])]
#     
#     
#         dates_crop = dates[date_crop:None if date_crop is None or date_crop==0 else -date_crop]
#         pairs_crop = pairs[pairs.ref.isin(dates_crop)&pairs.rep.isin(dates_crop)]
#     
#     #     refs = [ref for ref in pairs['ref'] if str(ref.date()) in dates]
#     #     reps = [rep for rep in pairs['rep'] if str(rep.date()) in dates]
#     
#         phase_turbo = xr.concat([(
#             (turbo.sel(date=ref).drop('date') if ref in turbo.date else fake) - \
#             (turbo.sel(date=rep).drop('date') if rep in turbo.date else fake)
#         ).assign_coords(pair=str(ref.date()) + ' ' + str(rep.date()), ref=ref, rep=rep) \
#         for ref, rep in zip(pairs_crop['ref'], pairs_crop['rep'])], dim='pair')
#     
#     #     phase_turbo = xr.concat([(turbo.sel(date=ref) - turbo.sel(date=rep))\
#     #                              .assign_coords(pair=str(ref.date()) + ' ' + str(rep.date()), ref=ref, rep=rep) \
#     #                              for ref, rep in \
#     #                              zip(pairs['ref'], pairs['rep'])],
#     #                        dim='pair')
#         #phase_turbo['ref'].values = pd.to_datetime(phase_turbo['ref'])
#         #phase_turbo['rep'].values = pd.to_datetime(phase_turbo['rep'])
#         return phase_turbo

#     def turbulence(self, phase, weight=None, date_crop=1, symmetrical=False):
#         import xarray as xr
#         import pandas as pd
# 
#         pairs, dates = self.get_pairs(phase, dates=True)
# 
#         turbos = []
#         for date in dates:
#             ref = pairs[pairs.ref==date]
#             rep = pairs[pairs.rep==date]
#             # calculate left and right pairs
#             count = min(len(ref), len(rep)) if symmetrical else None
#             #print (date, len(ref), len(rep), '=>', count)
#             if len(ref) < 1 or len(rep) < 1:
#                 # correction calculation is not possible
#                 continue
#             ref_data = phase.sel(pair=ref.pair.values).isel(pair=slice(None,count))
#             #print (ref_data)
#             rep_data = phase.sel(pair=rep.pair.values ).isel(pair=slice(None,count))
#             #print (rep_data)
#             if weight is not None:
#                 ref_weight = weight.sel(pair=ref_data.pair)
#                 rep_weight = weight.sel(pair=rep_data.pair)
#                 turbo = ((ref_data*ref_weight).mean('pair')/ref_weight.sum('pair') -\
#                          (rep_data*rep_weight).mean('pair')/rep_weight.sum('pair')) / 2
#                 del ref_weight, rep_weight
#             else:
#                 turbo = (ref_data.mean('pair') - rep_data.mean('pair')) / 2
#             del ref_data, rep_data
#             turbos.append(turbo.assign_coords(date=pd.to_datetime(date)))
#             del turbo
#         turbo = xr.concat(turbos, dim='date')
#         del turbos
# 
#         # empty grid
#         empty = xr.zeros_like(phase.isel(pair=0))
#         # convert dates to pairs
#         dates_crop = dates[date_crop:None if date_crop is None or date_crop==0 else -date_crop]
#         pairs_crop = pairs[pairs.ref.isin(dates_crop) & pairs.rep.isin(dates_crop)]  
#         phase_turbo = xr.concat([(
#             (turbo.sel(date=ref).drop('date') if ref in turbo.date else empty) - \
#             (turbo.sel(date=rep).drop('date') if rep in turbo.date else empty)
#         ).assign_coords(pair=str(ref.date()) + ' ' + str(rep.date()), ref=ref, rep=rep) \
#         for ref, rep in zip(pairs_crop['ref'], pairs_crop['rep'])], dim='pair')
#         del empty, dates_crop, pairs_crop
# 
#         return phase_turbo.rename('turbulence')

    def turbulence(self, phase, weight=None):
        import xarray as xr
        import pandas as pd

        pairs, dates = self.get_pairs(phase, dates=True)

        turbos = []
        for date in dates:
            ref = pairs[pairs.ref==date]
            rep = pairs[pairs.rep==date]
            #print (date, len(ref), len(rep))
            ref_data = phase.sel(pair=ref.pair.values)
            #print (ref_data)
            rep_data = phase.sel(pair=rep.pair.values)
            #print (rep_data)
            if weight is not None:
                ref_weight = weight.sel(pair=ref.pair.values)
                rep_weight = weight.sel(pair=rep.pair.values)
                turbo = xr.concat([ref_data*ref_weight, -rep_data*rep_weight], dim='pair').sum('pair')/\
                    xr.concat([ref_weight, rep_weight], dim='pair').sum('pair')
                del ref_weight, rep_weight
            else:
                turbo = xr.concat([ref_data, -rep_data], dim='pair').mean('pair')
            del ref_data, rep_data
            turbos.append(turbo.assign_coords(date=pd.to_datetime(date)))
            del turbo
        turbo = xr.concat(turbos, dim='date')
        del turbos

        phase_turbo = xr.concat([(turbo.sel(date=ref).drop('date') - turbo.sel(date=rep).drop('date'))\
                                 .assign_coords(pair=str(ref.date()) + ' ' + str(rep.date()), ref=ref, rep=rep) \
                          for ref, rep in zip(pairs['ref'], pairs['rep'])], dim='pair')

        return phase_turbo.rename('turbulence')
