# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_phasediff import Stack_phasediff

class Stack_multilooking(Stack_phasediff):

    #decimator = lambda da: da.coarsen({'y': 2, 'x': 2}, boundary='trim').mean()
    def decimator(self, resolution=60, grid=(1, 4), func='mean', debug=False):
        """
        Return function for pixel decimation to the specified output resolution.

        Parameters
        ----------
        resolution : int, optional
            DEM grid resolution in meters. The same grid is used for geocoded results output.
        grid : tuple, optional
            Grid size for pixel decimation in the format (vertical, horizontal).
        debug : bool, optional
            Boolean flag to print debug information.

        Returns
        -------
        callable
            Post-processing lambda function.

        Examples
        --------
        Decimate computed interferograms to default DEM resolution 60 meters:
        decimator = stack.decimator()
        stack.intf(pairs, func=decimator)
        """
        import numpy as np
        import dask
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')

        # special cases: scale factor should be 4*N or 2*N to prevent rounding issues
        # grid can be defined as [] or () or xarray dataarray
        grid_as_coeffs = isinstance(grid,(list, tuple))
        if   (grid_as_coeffs and grid==(1, 1)) or (not grid_as_coeffs and grid.x.diff('x')[0].item()) == 1:
            xscale0 = 4
        elif (grid_as_coeffs and grid==(1, 2)) or (not grid_as_coeffs and grid.x.diff('x')[0].item()) == 2:
            xscale0 = 2
        else:
            xscale0 = 1
        if debug:
            print (f'DEBUG: scale to square grid: xscale0={xscale0}')

        dy, dx = self.get_spacing(grid)
        yscale, xscale = int(np.round(resolution/dy)), int(np.round(resolution/dx/xscale0))
        if debug:
            print (f'DEBUG: ground pixel size in meters: y={dy:.1f}, x={dx:.1f}')
        if yscale <= 1 and xscale <= 1 and xscale0==1:
            # decimation impossible
            if debug:
                print (f'DEBUG: decimator = lambda da: da')
            return lambda da: da
        if debug:
            print (f"DEBUG: decimator = lambda da: da.coarsen({{'y': {yscale}, 'x': {xscale0*xscale}}}, boundary='trim').{func}()")

        # decimate function
        def decimator(da):
            import warnings
            # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
            warnings.filterwarnings('ignore')
            warnings.filterwarnings('ignore', module='dask')
            warnings.filterwarnings('ignore', module='dask.core')
            # workaround for Google Colab when we cannot save grids with x,y coordinate names
            # also supports geographic coordinates
            yname = [varname for varname in ['y', 'lat', 'a'] if varname in da.dims][0]
            xname = [varname for varname in ['x', 'lon', 'r'] if varname in da.dims][0]
            coarsen_args = {yname: yscale, xname: xscale0*xscale}
            #if debug:
            #    print (f"Decimate y variable '{yname}' for scale 1/{yscale} and x variable '{xname}' for scale 1/{xscale}")
            # avoid creating the large chunks
            with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                if func == 'mean':
                    return da.coarsen(coarsen_args, boundary='trim').mean()
                elif func == 'min':
                    return da.coarsen(coarsen_args, boundary='trim').min()
                elif func == 'max':
                    return da.coarsen(coarsen_args, boundary='trim').max()
                elif func == 'count':
                    return da.coarsen(coarsen_args, boundary='trim').count()
                elif func == 'sum':
                    return da.coarsen(coarsen_args, boundary='trim').sum()
                else:
                    raise ValueError(f"Unsupported function {func}. Should be 'mean','min','max','count', or 'sum'")

        # return callback function and set common chunk size
        return lambda da: decimator(da).chunk(self.chunksize)

    # coarsen = None disables downscaling and uses wavelength to filter
    # coarsen=1 disables downscaling and use coarsen/cutoff filter
    def multilooking(self, data, weight=None, wavelength=None, coarsen=None, debug=False):
        import xarray as xr
        import numpy as np
        import dask
        from dask_image.ndfilters import gaussian_filter as dask_gaussian_filter
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')
        # GMTSAR constant 5.3 defines half-gain at filter_wavelength
        # https://github.com/gmtsar/gmtsar/issues/411
        cutoff = 5.3

        # expand simplified definition
        if coarsen is not None and not isinstance(coarsen, (list, tuple, np.ndarray)):
            coarsen = (coarsen, coarsen)

        # allow this case to save the original grid resolution
        if wavelength is None and coarsen is None:
            return data

        # antialiasing (multi-looking) filter
        if wavelength is None:
            sigmas = [coarsen[0]/cutoff, coarsen[1]/cutoff]
            if debug:
                print (f'DEBUG: multilooking sigmas ({sigmas[0]:.2f}, {sigmas[1]:.2f}), for specified coarsen {coarsen}')
        else:
            dy, dx = self.get_spacing(data)
            #print ('DEBUG dy, dx', dy, dx)
            #sigmas = int(np.round(wavelength/dy/coarsen[0])), int(np.round(wavelength/dx))
            sigmas = [wavelength/cutoff/dy, wavelength/cutoff/dx]
            if debug:
                print (f'DEBUG: multilooking sigmas ({sigmas[0]:.2f}, {sigmas[1]:.2f}), for specified wavelength {wavelength:.1f}')

        # weighted and not weighted convolution on float and complex float data
        def apply_filter(data, weight, sigmas, truncate=2):
            import warnings
            # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
            warnings.filterwarnings('ignore')
            warnings.filterwarnings('ignore', module='dask')
            warnings.filterwarnings('ignore', module='dask.core')
            if np.issubdtype(data.dtype, np.complexfloating):
                parts = []
                for part in [data.real, data.imag]:
                    data_complex  = ((1j + part) * (weight if weight is not None else 1)).fillna(0)
                    conv_complex = dask_gaussian_filter(data_complex.data, sigmas, mode='reflect', truncate=truncate)
                    conv = conv_complex.real/conv_complex.imag
                    del data_complex, conv_complex
                    parts.append(conv)
                    del conv
                conv = parts[0] + 1j*parts[1]
                del parts
            else:
                # replace nan + 1j to to 0.+0.j
                data_complex  = ((1j + data) * (weight if weight is not None else 1)).fillna(0)
                conv_complex = dask_gaussian_filter(data_complex.data, sigmas, mode='reflect', truncate=truncate)
                conv = conv_complex.real/conv_complex.imag
                del data_complex, conv_complex
            return conv

        if isinstance(data, xr.Dataset):
            dims = data[list(data.data_vars)[0]].dims
        else:
            dims = data.dims

        if len(dims) == 2:
            stackvar = None
        else:
            stackvar = dims[0]

        if weight is not None and len(data.dims) == len(weight.dims):
            # single 2D grid processing
            assert data.shape == weight.shape, 'ERROR: multilooking data and weight variables have different shape'
        elif weight is not None and len(data.dims) + 1 == len(weight.dims):
            # stack of 2D grids processing
            assert data.shape[1:] == weight.shape, 'ERROR: multilooking data and weight variables have different shape'

        stack =[]
        for ind in range(len(data[stackvar]) if stackvar is not None else 1):
            if isinstance(data, xr.Dataset):
                data_convs = []
                for key in data.data_vars:
                    conv = apply_filter(data[key][ind] if stackvar is not None else data[key],
                                        weight[ind]    if stackvar is not None and weight is not None else weight,
                                        sigmas)
                    data_conv = xr.DataArray(conv, dims=data[key].dims[1:] if stackvar is not None else data[key].dims, name=data[key].name)
                    del conv
                    data_convs.append(data_conv)
                    del data_conv
                stack.append(xr.merge(data_convs))
                del data_convs
            else:
                conv = apply_filter(data[ind]   if stackvar is not None else data,
                                    weight[ind] if stackvar is not None and weight is not None else weight,
                                    sigmas)
                data_conv = xr.DataArray(conv, dims=data.dims[1:] if stackvar is not None else data.dims, name=data.name)
                del conv
                stack.append(data_conv)
                del data_conv

        if stackvar is not None:
            #print ('3D')
            ds = xr.concat(stack, dim=stackvar).assign_coords(data.coords)
        else:
            #print ('2D')
            ds = stack[0].assign_coords(data.coords)
        del stack

        # it works faster when we prevent small chunks usage
        if stackvar is not None:
            chunksizes = (1, self.chunksize, self.chunksize)
        else:
            chunksizes = (self.chunksize, self.chunksize)
        if coarsen is not None:
            # coarse grid to square cells
            return ds.coarsen({'y': coarsen[0], 'x': coarsen[1]}, boundary='trim').mean().chunk(chunksizes)
        return ds.chunk(chunksizes)
