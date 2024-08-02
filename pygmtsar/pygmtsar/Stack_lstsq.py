# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2024, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_tidal import Stack_tidal
from .tqdm_dask import tqdm_dask

class Stack_lstsq(Stack_tidal):

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
        pairs : pandas.DataFrame or Xarray object
            DataFrame containing interferogram date pairs.

        Returns
        -------
        numpy.ndarray
            Matrix with one row for every interferogram and one column for every date.
            Each element in the matrix is an integer, with 1 representing that the date
            is between the corresponding interferogram's reference and repeat dates, and
            0 otherwise.
        """
        return (self.get_pairs_matrix(pairs)>=0).astype(int)

    def lstsq_matrix_edge(self, pairs):
        """
        Create an edge matrix for use in the least squares computation based on interferogram date pairs.

        Parameters
        ----------
        pairs : pandas.DataFrame or Xarray object
            DataFrame containing interferogram date pairs.

        Returns
        -------
        numpy.ndarray
            Matrix with one row for every interferogram and one column for every date.
            Each element in the matrix is an integer, with 1 representing the end date,
            -1 the start date and 0 otherwise.
        """
        print ('NOTE: this function is not used in the code and created for test purposes only.')
        import numpy as np
        return np.nan_to_num(self.get_pairs_matrix(pairs)).astype(int)

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

        if not 'stack' in data.dims:
            chunks_z, chunks_y, chunks_x = data.chunks if data.chunks is not None else np.inf, np.inf, np.inf
            if np.max(chunks_y) > self.netcdf_chunksize or np.max(chunks_x) > self.netcdf_chunksize:
                print (f'Note: data chunk size ({np.max(chunks_y)}, {np.max(chunks_x)}) is too large for stack processing')
                chunks_y = chunks_x = self.netcdf_chunksize//2
                print (f'Note: auto tune data chunk size to a half of NetCDF chunk: ({chunks_y}, {chunks_x})')
            chunks = {'y': chunks_y, 'x': chunks_x}
        else:
            chunks_z, chunks_stack = data.chunks if data.chunks is not None else np.inf, np.inf
            if np.max(chunks_stack) > self.chunksize1d:
                print (f'Note: data chunk size ({np.max(chunks_stack)} is too large for stack processing')
                # 1D chunk size can be defined straightforward
                chunks_stack = self.chunksize1d
                print (f'Note: auto tune data chunk size to 1D chunk: ({chunks_stack})')
            chunks = {'stack': chunks_stack}
        data = data.chunk(chunks)

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
            #if np.max(chunks_y) > self.netcdf_chunksize or np.max(chunks_x) > self.netcdf_chunksize:
            #    print ('Note: auto tune weight chunk size to a half of NetCDF chunk')
            #    weight = weight.chunk({'y': chunks_y, 'x': chunks_x})
            weight = weight.chunk({'y': chunks_y, 'x': chunks_x})
        elif 'stack' in weight.dims:
            # this case should be processed inside lstq_block function
            assert weight.shape == data.shape, 'ERROR: data and weight dataarrays should have the same dimensions'
            weight = weight.chunk({'stack': chunks_stack})
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

        def lstq_block(ys, xs, stacks=None):
            if stacks is None:
                # 3D array
                data_block = data.isel(y=ys, x=xs).compute(n_workers=1).values.transpose(1,2,0)
            else:
                # 2D array
                data_block = data.isel(stack=stacks).compute(n_workers=1).values.transpose(1,0)
            if np.isnan(data_block).all():
                # do not process an empty block
                if stacks is None:
                    shape = (matrix.shape[1], *data_block.shape[:2])
                else:
                    shape = (matrix.shape[1], data_block.shape[0])
                return np.full(shape, np.nan, dtype=np.float32)
            # weights can be defined by multiple ways or be set to None
            if isinstance(weight, xr.DataArray):
                if stacks is None:
                    # 3D array
                    weight_block = weight.isel(y=ys, x=xs).compute(n_workers=1).values.transpose(1,2,0)
                else:
                    weight_block = weight.isel(stack=stacks).compute(n_workers=1).values.transpose(1,0)
                # weight=1 is not allowed for the used weighted least squares calculation function 
                weight_block = np.where(weight_block>=1, 1, weight_block)
                # Vectorize vec_lstsq
                vec_lstsq = np.vectorize(lambda x, w: self.lstsq1d(x, w, matrix, cumsum), signature='(n),(n)->(m)')
                # Apply vec_lstsq to data_block and weight_block and revert the original dimensions order
                if stacks is None:
                    block = vec_lstsq(data_block, weight_block).transpose(2,0,1)
                else:
                    block = vec_lstsq(data_block, weight_block).transpose(1,0)
                del weight_block, vec_lstsq
            else:
                # Vectorize vec_lstsq
                vec_lstsq = np.vectorize(lambda x: self.lstsq1d(x, weight, matrix, cumsum), signature='(n)->(m)')
                # Apply vec_lstsq to data_block and weight_block and revert the original dimensions order
                if stacks is None:
                    block = vec_lstsq(data_block).transpose(2,0,1)
                else:
                    block = vec_lstsq(data_block).transpose(1,0)
                del vec_lstsq
            del data_block
            return block

        # split to chunks
        # use indices instead of the coordinate values to prevent the weird error raising occasionally:
        # "Reindexing only valid with uniquely valued Index objects"
        if not 'stack' in data.dims:
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
        else:
            chunks_z, chunks_stack = data.chunks
            stacks_blocks = np.array_split(np.arange(data['stack'].size), np.cumsum(chunks_stack)[:-1])
            blocks_total = []
            for stacks_block in stacks_blocks:
                block = dask.array.from_delayed(dask.delayed(lstq_block)(None, None, stacks_block),
                                                shape=(len(dates), stacks_block.size),
                                                dtype=np.float32)
                blocks_total.append(block)
                del block
            model = dask.array.block(blocks_total)
            del blocks_total
            coords = {'date': pd.to_datetime(dates), 'stack': data['stack']}
            
        model = xr.DataArray(model, coords=coords).rename('displacement')

        return model

    def rmse(self, data, solution, weight=None):
        """
        Calculate difference between pairs and dates
    
        Use to calculate solution vs pair unwrapped phases difference as
        diff = sbas.stack_vs_cube(phase_unwrap, solution)
        error = np.sqrt(sbas.wrap(diff)**2).sum('pair')
        """
        import numpy as np
        import xarray as xr
        import pandas as pd

        multi_index = None
        if 'stack' in data.dims and isinstance(data.coords['stack'].to_index(), pd.MultiIndex):
            assert 'stack' in solution.dims, 'ERROR: "solution" must be stacked consistently with "data".'
            if 'y' in data.coords and 'x' in data.coords:
                multi_index_names = ['y', 'x']
            elif 'lat' in data.coords and 'lon' in data.coords:
                multi_index_names = ['lat', 'lon']
            multi_index = pd.MultiIndex.from_arrays([data.y.values, data.x.values], names=multi_index_names)
            data = data.reset_index('stack')
            solution = solution.reset_index('stack')
            if weight is not None:
                assert 'stack' in weight.dims, 'ERROR: "weight" must be stacked consistently with "data".'
                weight = weight.reset_index('stack')
        elif not 'stack' in data.dims:
            assert not 'stack' in solution.dims, 'ERROR: "solution" must be stacked consistently with "data".'
            if weight is not None:
                assert not 'stack' in weight.dims, 'ERROR: "weight" must be stacked consistently with "data".'
        
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
        if weight is not None:
            return np.sqrt((weight * error).sum('pair') / weight.sum('pair') / len(pairs)).rename('rmse')
        rmse = np.sqrt((error).sum('pair') / len(pairs)).rename('rmse')
        if multi_index is not None:
            return rmse.assign_coords(stack=multi_index)
        return rmse

    def plot_displacement(self, data, caption='Cumulative LOS Displacement, [rad]',
                          quantile=None, vmin=None, vmax=None, symmetrical=False, aspect=None, **kwargs):
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt

        if 'stack' in data.dims and isinstance(data.coords['stack'].to_index(), pd.MultiIndex):
            data = data.unstack('stack')

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
                           quantile=None, vmin=None, vmax=None, symmetrical=False, **kwargs):
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt

        if 'stack' in data.dims and isinstance(data.coords['stack'].to_index(), pd.MultiIndex):
            data = data.unstack('stack')

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
        if self.is_ra(data):
            fg.set_axis_labels('Range', 'Azimuth')
        fg.set_ticks(max_xticks=nbins, max_yticks=nbins)
        fg.fig.suptitle(caption, y=y)
        
        self.plots_AOI(fg, **kwargs)
        self.plots_POI(fg, **kwargs)

    def plot_displacements_los_mm(self, data, caption='Cumulative LOS Displacement, [mm]', cols=4, size=4, nbins=5, aspect=1.2, y=1.05,
                           quantile=None, vmin=None, vmax=None, symmetrical=False, **kwargs):
        self.plot_displacements(self.los_displacement_mm(data),
                                caption=caption, cols=cols, size=size, nbins=nbins, aspect=aspect, y=y,
                                quantile=quantile, vmin=vmin, vmax=vmax, symmetrical=symmetrical, **kwargs)

    def plot_rmse(self, data, caption='RMSE, [rad]', cmap='turbo',
                  quantile=None, vmin=None, vmax=None, symmetrical=False, **kwargs):
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')

        if 'stack' in data.dims and isinstance(data.coords['stack'].to_index(), pd.MultiIndex):
            data = data.unstack('stack')

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
        data.plot.imshow(cmap=cmap, vmin=vmin, vmax=vmax)
        self.plot_AOI(**kwargs)
        self.plot_POI(**kwargs)
        plt.title(caption)

    def plot_rmse_los_mm(self, data, caption='RMSE, [mm]', cmap='turbo',
                  quantile=None, vmin=None, vmax=None, symmetrical=False, **kwargs):
        self.plot_rmse(abs(self.los_displacement_mm(data)),
                       caption=caption, cmap=cmap,
                       quantile=quantile, vmin=vmin, vmax=vmax, symmetrical=symmetrical, **kwargs)
