# ----------------------------------------------------------------------------
# InSAR.dev
# 
# This file is part of the InSAR.dev project: https://InSAR.dev
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_phasediff import Stack_phasediff
from .utils import utils

class Stack_multilooking(Stack_phasediff):

    #decimator = lambda da: da.coarsen({'y': 2, 'x': 2}, boundary='trim').mean()
    def decimator(self, grid, resolution=60, func='mean', debug=False):
        """
        Return function for pixel decimation to the specified output resolution.

        Parameters
        ----------
        grid : xarray object
            Grid to define the spacing.
        resolution : int, optional
            DEM grid resolution in meters. The same grid is used for geocoded results output.
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

        dy, dx = self.get_spacing(grid)
        yscale, xscale = int(np.round(resolution/dy)), int(np.round(resolution/dx))
        if debug:
            print (f'DEBUG: ground pixel size in meters: y={dy:.1f}, x={dx:.1f}')
        if yscale <= 1 and xscale <= 1:
            # decimation impossible
            if debug:
                print (f'DEBUG: decimator = lambda da: da')
            return lambda da: da
        if debug:
            print (f"DEBUG: decimator = lambda da: da.coarsen({{'y': {yscale}, 'x': {xscale}}}, boundary='trim').{func}()")

        # decimate function
        def decimator(da):
            import warnings
            # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
            warnings.filterwarnings('ignore')
            warnings.filterwarnings('ignore', module='dask')
            warnings.filterwarnings('ignore', module='dask.core')
            # unstack data if needed
            if 'stack' in da.dims:
                # .unstack() is too slow on lazy grids in some of Xarray/Dask versions
                da = da.compute().unstack('stack')
            # workaround for Google Colab when we cannot save grids with x,y coordinate names
            # also supports geographic coordinates
            yname = [varname for varname in ['y', 'lat', 'a'] if varname in da.dims][0]
            xname = [varname for varname in ['x', 'lon', 'r'] if varname in da.dims][0]
            coarsen_args = {yname: yscale, xname: xscale}
            # calculate coordinate offsets to align coarsened grids
            y0 = self.calculate_coarsen_start(da, yname, yscale)
            x0 = self.calculate_coarsen_start(da, xname, xscale)
            # avoid creating the large chunks
            with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                #if func not in ['mean', 'min', 'max', 'count', 'sum']:
                #    raise ValueError(f"Unsupported function {func}. Should be 'mean','min','max','count', or 'sum'")
                # return getattr(da.coarsen(coarsen_args, boundary='trim'), func)()\
                #        .chunk({yname: self.chunksize, xname: self.chunksize})
                return getattr(da.isel({yname: slice(y0, None), xname: slice(x0, None)})\
                       .coarsen(coarsen_args, boundary='trim'), func)()\
                       .chunk({yname: self.chunksize, xname: self.chunksize})

        # return callback function and set common chunk size
        return lambda da: decimator(da)

    def multilooking(self, data, weight=None, wavelength=None, coarsen=None, debug=False):
        import xarray as xr
        import numpy as np
        import dask
    
        # GMTSAR constant 5.3 defines half-gain at filter_wavelength
        cutoff = 5.3
    
        # Expand simplified definition of coarsen
        coarsen = (coarsen, coarsen) if coarsen is not None and not isinstance(coarsen, (list, tuple, np.ndarray)) else coarsen
    
        # no-op, processing is needed
        if wavelength is None and coarsen is None:
            return data
    
        # calculate sigmas based on wavelength or coarsen
        if wavelength is not None:
            dy, dx = self.get_spacing(data)
            sigmas = [wavelength / cutoff / dy, wavelength / cutoff / dx]
            if debug:
                print(f'DEBUG: multilooking sigmas ({sigmas[0]:.2f}, {sigmas[1]:.2f}), wavelength {wavelength:.1f}')
        else:
            sigmas = [coarsen[0] / cutoff, coarsen[1] / cutoff]
            if debug:
                print(f'DEBUG: multilooking sigmas ({sigmas[0]:.2f}, {sigmas[1]:.2f}), coarsen {coarsen}')

        if isinstance(data, xr.Dataset):
            dims = data[list(data.data_vars)[0]].dims
        else:
            dims = data.dims

        if len(dims) == 2:
            stackvar = None
        else:
            stackvar = dims[0]
        #print ('stackvar', stackvar)

        if weight is not None:
            # for InSAR processing expect 2D weights
            assert isinstance(weight, xr.DataArray) and len(weight.dims)==2, \
                'ERROR: multilooking weight should be 2D DataArray'
        
        if weight is not None and len(data.dims) == len(weight.dims):
            #print ('2D check shape weighted')
            # single 2D grid processing
            if isinstance(data, xr.Dataset):
                for varname in data.data_vars:
                    assert data[varname].shape == weight.shape, \
                        f'ERROR: multilooking data[{varname}] and weight variables have different shape'
            else:
                assert data.shape == weight.shape, 'ERROR: multilooking data and weight variables have different shape'
        elif weight is not None and len(data.dims) == len(weight.dims) + 1:
            #print ('3D check shape weighted')
            # stack of 2D grids processing
            if isinstance(data, xr.Dataset):
                for varname in data.data_vars:
                    assert data[varname].shape[1:] == weight.shape, \
                        f'ERROR: multilooking data[{varname}] slice and weight variables have different shape \
                        ({data[varname].shape[1:]} vs {weight.shape})'
            else:
                assert data.shape[1:] == weight.shape, f'ERROR: multilooking data slice and weight variables have different shape \
                ({data.shape[1:]} vs {weight.shape})'

        # process a slice of dataarray
        def process_slice(slice_data):
            conv = self.nanconvolve2d_gaussian(slice_data, weight, sigmas)
            return xr.DataArray(conv, dims=slice_data.dims, name=slice_data.name)

        # process stack of dataarray slices
        def process_slice_var(dataarray):    
            if stackvar:
                stack = [process_slice(dataarray[ind]) for ind in range(len(dataarray[stackvar]))]
                return xr.concat(stack, dim=stackvar).assign_coords(dataarray.coords)
            else:
                return process_slice(dataarray).assign_coords(dataarray.coords)

        if isinstance(data, xr.Dataset):
            ds = xr.Dataset({varname: process_slice_var(data[varname]) for varname in data.data_vars})
        else:
            ds = process_slice_var(data)
    
        # Set chunk size
        chunksizes = {'y': self.chunksize, 'x': self.chunksize}

        if coarsen:
            # calculate coordinate offsets to align coarsened grids
            y0 = self.calculate_coarsen_start(ds, 'y', coarsen[0])
            x0 = self.calculate_coarsen_start(ds, 'x', coarsen[1])
            ds = ds.isel({'y': slice(y0, None), 'x': slice(x0, None)})\
                     .coarsen({'y': coarsen[0], 'x': coarsen[1]}, boundary='trim')\
                     .mean()

        return self.spatial_ref(ds.chunk(chunksizes), data)
