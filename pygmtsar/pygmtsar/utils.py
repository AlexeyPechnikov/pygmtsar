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
                *(da.broadcast_arrays(data, weight) if weight is not None else [data]),
                depth=depth,
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
