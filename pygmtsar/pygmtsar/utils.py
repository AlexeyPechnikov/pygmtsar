# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
class utils():

#     @staticmethod
#     def regression(phase, topo, fit_intercept=True):
#         import numpy as np
#         import xarray as xr
#         from sklearn.linear_model import LinearRegression
#         from sklearn.pipeline import make_pipeline
#         from sklearn.preprocessing import StandardScaler
# 
#         # define on the same grid
#         topo = topo.reindex_like(phase, method='nearest')
#     
#         # build prediction model with or without plane removal (fit_intercept)
#         regr = make_pipeline(StandardScaler(), LinearRegression(fit_intercept=fit_intercept))
# 
#         topo_values = topo.values.ravel()
#         phase_values = phase.values.ravel()
#         nanmask = np.isnan(a) | np.isnan(b)
# 
#         # fit on non-NaN values only and predict on the full grid
#         phase_topo = regr.fit(np.column_stack([topo_values[~nanmask]]),
#                               np.column_stack([phase_values[~nanmask]]))\
#                     .predict(np.column_stack([topo_values])).reshape(phase.shape)
#         return xr.DataArray(phase_topo, coords=phase.coords)

    # Xarray's interpolation can be inefficient for large grids;
    # this custom function handles the task more effectively.
    @staticmethod
    def interp2d_like(data, grid, method='cubic', **kwargs):
        """
        Efficiently interpolate a 2D array using OpenCV interpolation methods.
        
        Args:
            data (xarray.DataArray): The input data array.
            grid (xarray.DataArray): The grid to interpolate onto.
            method (str): Interpolation method ('nearest', 'linear', 'cubic' or 'lanczos').
            **kwargs: Additional arguments for interpolation.
    
        Returns:
            xarray.DataArray: The interpolated data.
        """
        import cv2
        import numpy as np
        import xarray as xr
        import dask.array as da
        dims = grid.dims[-2:]
        dim1, dim2 = dims
        coords = {dim1: grid[dim1], dim2: grid[dim2]}
        #print ('coords', coords)
    
        # Define interpolation method
        if method == 'nearest':
            interpolation = cv2.INTER_NEAREST
        elif method == 'linear':
            interpolation = cv2.INTER_LINEAR
        elif method == 'cubic':
            interpolation = cv2.INTER_CUBIC
        elif method == 'lanczos':
            interpolation = cv2.INTER_LANCZOS4
        else:
            raise ValueError(f"Unsupported interpolation {method}. Should be 'nearest', 'linear', 'cubic' or 'lanczos'")

        # TBD: can be added to the function parameters
        borderMode = cv2.BORDER_REFLECT
    
        # define interpolation function using outer variable data
        def interpolate_chunk(out_chunk1, out_chunk2, dim1, dim2, interpolation, borderMode, **kwargs):
            d1 = float(data[dim1].diff(dim1)[0])
            d2 = float(data[dim2].diff(dim2)[0])
    
            # select the chunk from data with some padding
            chunk = data.sel({
                dim1: slice(out_chunk1[0] - 3 * d1, out_chunk1[-1] + 3 * d1),
                dim2: slice(out_chunk2[0] - 3 * d2, out_chunk2[-1] + 3 * d2)
            }).compute(n_workers=1)
    
            # Create grid for interpolation
            dst_grid_x, dst_grid_y = np.meshgrid(out_chunk2, out_chunk1)
    
            # map destination grid coordinates to source pixel indices
            src_x_coords = np.interp(
                dst_grid_x.ravel(),
                chunk[dim2].values,
                np.arange(len(chunk[dim2]))
            )
            src_y_coords = np.interp(
                dst_grid_y.ravel(),
                chunk[dim1].values,
                np.arange(len(chunk[dim1]))
            )
    
            # reshape the coordinates for remap
            src_x_coords = src_x_coords.reshape(dst_grid_x.shape).astype(np.float32)
            src_y_coords = src_y_coords.reshape(dst_grid_y.shape).astype(np.float32)
    
            # interpolate using OpenCV
            dst_grid = cv2.remap(
                chunk.values.astype(np.float32),
                src_x_coords,
                src_y_coords,
                interpolation=interpolation,
                borderMode=borderMode
            )
            return dst_grid
    
        # define chunk sizes
        chunk_sizes = grid.chunks[-2:] if hasattr(grid, 'chunks') else (data.sizes[dim1], data.sizes[dim2])
    
        # create dask array for parallel processing
        grid_y = da.from_array(grid[dim1].values, chunks=chunk_sizes[0])
        grid_x = da.from_array(grid[dim2].values, chunks=chunk_sizes[1])
    
        # Perform interpolation
        dask_out = da.blockwise(
            interpolate_chunk,
            'yx',
            grid_y, 'y',
            grid_x, 'x',
            dtype=data.dtype,
            dim1=dim1,
            dim2=dim2,
            interpolation=interpolation,
            borderMode=borderMode,
            **kwargs
        )
    
        da_out = xr.DataArray(dask_out, coords=coords, dims=dims).rename(data.name)
    
        # Append all the input coordinates
        return da_out.assign_coords({k: v for k, v in data.coords.items() if k not in coords})

#     # Xarray's interpolation can be inefficient for large grids;
#     # this custom function handles the task more effectively.
#     @staticmethod
#     def interp2d_like(data, grid, method='cubic', **kwargs):
#         import xarray as xr
#         import dask.array as da
#         import os
#         import warnings
#         # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
#         warnings.filterwarnings('ignore')
#         warnings.filterwarnings('ignore', module='dask')
#         warnings.filterwarnings('ignore', module='dask.core')
# 
#         # detect dimensions and coordinates for 2D or 3D grid
#         dims = grid.dims[-2:]
#         dim1, dim2 = dims
#         coords = {dim1: grid[dim1], dim2: grid[dim2]}
#         #print (f'dims: {dims}, coords: {coords}')
# 
#         # use outer variable data
#         def interpolate_chunk(out_chunk1, out_chunk2, dim1, dim2, method, **kwargs):
#             d1, d2 = float(data[dim1].diff(dim1)[0]), float(data[dim2].diff(dim2)[0])
#             #print ('d1, d2', d1, d2)
#             chunk = data.sel({
#                                 dim1: slice(out_chunk1[0]-2*d1, out_chunk1[-1]+2*d1),
#                                 dim2: slice(out_chunk2[0]-2*d2, out_chunk2[-1]+2*d2)
#                                 }).compute(n_workers=1)
#             #print ('chunk', chunk)
#             out = chunk.interp({dim1: out_chunk1, dim2: out_chunk2}, method=method, **kwargs)
#             del chunk
#             return out
# 
#         chunk_sizes = grid.chunks[-2:] if hasattr(grid, 'chunks') else (self.chunksize, self.chunksize)
#         # coordinates are numpy arrays
#         grid_y = da.from_array(grid[dim1].values, chunks=chunk_sizes[0])
#         grid_x = da.from_array(grid[dim2].values, chunks=chunk_sizes[1])
# 
#         dask_out = da.blockwise(
#             interpolate_chunk,
#             'yx',
#             grid_y, 'y',
#             grid_x, 'x',
#             dtype=data.dtype,
#             dim1=dim1,
#             dim2=dim2,
#             method=method,
#             **kwargs
#         )
#         da_out = xr.DataArray(dask_out, coords=coords, dims=dims).rename(data.name)
#         del dask_out
#         # append all the input coordinates
#         return da_out.assign_coords({k: v for k, v in data.coords.items() if k not in coords})

    @staticmethod
    def nanconvolve2d_gaussian(data,
                        weight=None,
                        sigma=None,
                        mode='reflect',
                        truncate=4.0):
        import numpy as np
        import xarray as xr
    
        if sigma is None:
            return data
    
        if not isinstance(sigma, (list, tuple, np.ndarray)):
            sigma = (sigma, sigma)
        depth = [np.ceil(_sigma * truncate).astype(int) for _sigma in sigma]
        #print ('sigma', sigma, 'depth', depth)
    
        # weighted Gaussian filtering for real floats with NaNs
        def nanconvolve2d_gaussian_floating_dask_chunk(data, weight=None, **kwargs):
            import numpy as np
            from scipy.ndimage import gaussian_filter
            assert not np.issubdtype(data.dtype, np.complexfloating)
            assert np.issubdtype(data.dtype, np.floating)
            if weight is not None:
                assert not np.issubdtype(weight.dtype, np.complexfloating)
                assert np.issubdtype(weight.dtype, np.floating)
            # replace nan + 1j to to 0.+0.j
            data_complex  = (1j + data) * (weight if weight is not None else 1)
            conv_complex = gaussian_filter(np.nan_to_num(data_complex, 0), **kwargs)
            #conv = conv_complex.real/conv_complex.imag
            # to prevent "RuntimeWarning: invalid value encountered in divide" even when warning filter is defined
            conv = np.where(conv_complex.imag == 0, np.nan, conv_complex.real/(conv_complex.imag + 1e-17))
            del data_complex, conv_complex
            return conv
    
        def nanconvolve2d_gaussian_dask_chunk(data, weight=None, **kwargs):
            import numpy as np
            if np.issubdtype(data.dtype, np.complexfloating):
                #print ('complexfloating')
                real = nanconvolve2d_gaussian_floating_dask_chunk(data.real, weight, **kwargs)
                imag = nanconvolve2d_gaussian_floating_dask_chunk(data.imag, weight, **kwargs)
                conv = real + 1j*imag
                del real, imag
            else:
                #print ('floating')
                conv = nanconvolve2d_gaussian_floating_dask_chunk(data.real, weight, **kwargs)
            return conv
    
        # weighted Gaussian filtering for real or complex floats
        def nanconvolve2d_gaussian_dask(data, weight, **kwargs):
            import dask.array as da
            # ensure both dask arrays have the same chunk structure
            # use map_overlap with the custom function to handle both arrays
            return da.map_overlap(
                nanconvolve2d_gaussian_dask_chunk,
                *([data, weight] if weight is not None else [data]),
                depth={0: depth[0], 1: depth[1]},
                boundary='none',
                dtype=data.dtype,
                meta=data._meta,
                **kwargs
            )
            
        return xr.DataArray(nanconvolve2d_gaussian_dask(data.data,
                                     weight.data if weight is not None else None,
                                     sigma=sigma,
                                     mode=mode,
                                     truncate=truncate),
                            coords=data.coords,
                            name=data.name)

    @staticmethod
    def histogram(data, bins, range):
        """
        hist, bins = utils.histogram(corr60m.mean('pair'), 10, (0,1))
        print ('bins', bins)
        hist = dask.compute(hist)[0]
        print (hist)
        """
        import dask
        return dask.array.histogram(data, bins=bins, range=range)

    @staticmethod
    def corrcoef(data, mask_diagonal=False):
        """
        Calculate the correlation coefficients matrix for a stack of interferograms.

        Parameters:
        - data (xarray.DataArray): A 3D input grid representing a stack of interferograms.
                                  The dimensions should include time, 'y', and 'x'.
        - mask_diagonal (bool): If True, the diagonal elements of the correlation matrix
                                will be set to NaN to exclude self-correlation.

        Returns:
        xarray.DataArray: The cross-correlation matrix.
        
        Example:
        # plot variance-covariance matrix
        corr_stack = corr60m.mean('pair')
        corr = utils.corrcoef(unwrap.phase.where(corr_stack.where(corr_stack>0.7)))
        corr.where((corr)>0.7).plot.imshow(vmin=-1, vmax=1)
        plt.xticks(rotation=90)
        plt.show()
        
        # find the most correlated interferograms to exclude
        corr = utils.corrcoef(unwrap.phase.where(corr_stack.where(corr_stack>0.7)), mask_diagonal=True)
        df = corr.where(corr>0.7).to_dataframe('').reset_index().dropna()
        exclude = np.unique(np.concatenate([df.ref, df.rep]))
        exclude
        """
        import xarray as xr
        import dask
        import numpy as np
        assert len(data.dims) == 3, 'ERROR: expected 3D input grid'
        stackdim = data.dims[0]
        stackvals = data[stackdim].values
        corr = dask.array.corrcoef(data.stack(points=['y', 'x']).dropna(dim='points').data).round(2)
        corr = xr.DataArray((corr), coords={'ref': stackvals, 'rep': stackvals}).rename('corr')
        corr['ref'] = corr['ref'].astype(str)
        corr['rep'] = corr['rep'].astype(str)
        if mask_diagonal:
            # set diagonal elements to NaN
            diagonal_mask = np.eye(len(stackvals), dtype=bool)
            return corr.where(~diagonal_mask, np.nan)
        return corr

    @staticmethod
    def binary_erosion(data, *args, **kwargs):
        """
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_erosion.html
        """
        import xarray as xr
        from scipy.ndimage import binary_erosion
        array = binary_erosion(data.values, *args, **kwargs)
        return xr.DataArray(array, coords=data.coords, dims=data.dims, attrs=data.attrs)

    @staticmethod
    def binary_dilation(data, *args, **kwargs):
        """
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_dilation.html
        """
        import xarray as xr
        from scipy.ndimage import binary_dilation
        array = binary_dilation(data.values, *args, **kwargs)
        return xr.DataArray(array, coords=data.coords, dims=data.dims, attrs=data.attrs)

    @staticmethod
    def binary_opening(data, *args, **kwargs):
        """
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_opening.html
        
        corrmask = utils.binary_closing(corrmask, structure=np.ones((10,10)))
        corrmask = utils.binary_opening(corrmask, structure=np.ones((10,10)))
        """
        import xarray as xr
        from scipy.ndimage import binary_opening
        array = binary_opening(data.values, *args, **kwargs)
        return xr.DataArray(array, coords=data.coords, dims=data.dims, attrs=data.attrs)

    @staticmethod
    def binary_closing(data, *args, **kwargs):
        """
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_opening.html
        
        corrmask = utils.binary_closing(corrmask, structure=np.ones((10,10)))
        corrmask = utils.binary_opening(corrmask, structure=np.ones((10,10)))
        """
        import xarray as xr
        from scipy.ndimage import binary_closing
        array = binary_closing(data.values, *args, **kwargs)
        return xr.DataArray(array, coords=data.coords, dims=data.dims, attrs=data.attrs)
