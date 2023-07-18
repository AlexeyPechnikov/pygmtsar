# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2022, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------

class datagrid:
    """
    A class representing a data grid.

    Attributes
    ----------
    chunksize : int
        The chunk size for data compression. Default is 512.
    engine : str
        The engine used for NetCDF file operations. Default is 'h5netcdf'.
    complevel : int
        The compression level for data compression. Default is 3.
    noindex : np.uint32
        The NODATA index value for transform matrices.
    """
    import numpy as np

    # NetCDF options, see https://docs.xarray.dev/en/stable/user-guide/io.html#zarr-compressors-and-filters
    chunksize = 1024
    engine = 'h5netcdf'
    complevel = 3
    # NODATA index value for transform matrices
    # TODO: use the same datatype for the matrices to allow 64 bit datatype
    noindex = np.uint32(-1)

    # define lost class variables due to joblib via arguments
    def compression(self, shape=None, complevel=None, chunksize=None):
        """
        Return the compression options for a data grid.

        Parameters
        ----------
        shape : tuple, list, np.ndarray, optional
            The shape of the data grid. Required if chunksize is less than grid dimension sizes. Default is None.
        complevel : int, optional
            The compression level for data compression. If not specified, the class attribute complevel is used.
        chunksize : int or tuple, optional
            The chunk size for data compression. If not specified, the class attribute chunksize is used.

        Returns
        -------
        dict
            A dictionary containing the compression options for the data grid.

        Examples
        --------
        Get the compression options for a data grid with shape (1000, 1000):

        >>> compression(shape=(1000, 1000))
        {'zlib': True, 'complevel': 3, 'chunksizes': (512, 512)}

        Get the compression options for a data grid with chunksize 256:

        >>> compression(chunksize=256)
        {'zlib': True, 'complevel': 3, 'chunksizes': (256, 256)}
        """
        import numpy as np

        if complevel is None:
            complevel = self.complevel
        if chunksize is None:
            chunksize = self.chunksize
        assert chunksize is not None, 'compression() chunksize is None'
        if isinstance(chunksize, (tuple, list, np.ndarray)):
            # use as is, it can be 2D or 3D grid (even 1D while it is not used for now)
            if shape is not None:
                assert len(shape) == len(chunksize), 'ERROR: defined shape and chunksize dimensions are not equal'
                chunksizes = tuple([chunksize[dim] if chunksize[dim]<shape[dim] else shape[dim] for dim in range(len(shape))])
            else:
                chunksizes = chunksize
        else:
            # suppose 2D grid
            if shape is not None:
                chunksizes=(chunksize if chunksize<shape[0] else shape[0], chunksize if chunksize<shape[1] else shape[1])
            else:
                chunksizes=(chunksize, chunksize)
        return dict(zlib=True, complevel=complevel, chunksizes=chunksizes)

    @staticmethod
    def is_ra(grid):
        """
        Checks if the given grid is in radar coordinates.

        Parameters
        ----------
        grid : xarray.DataArray
            The grid to check.

        Returns
        -------
        bool
            True if the grid is in radar coordinates, False otherwise.
        """
        dims = grid.dims
        if 'y' in dims and 'x' in dims:
            return True
        return False

    @staticmethod
    def is_geo(grid):
        """
        Checks if the given grid is in geographic coordinates.

        Parameters
        ----------
        grid : xarray.DataArray
            The grid to check.

        Returns
        -------
        bool
            True if the grid is in geographic coordinates, False otherwise.
        """
        dims = grid.dims
        if 'lat' in dims and 'lon' in dims:
            return True
        return False

    def as_geo(self, da):
        """
        Add geospatial attributes (CRS and spatial dimensions) to allow raster operations using RioXarray.

        Parameters
        ----------
        da : xarray.DataArray
            The input 2D or 3D grid to be converted to geospatial.

        Returns
        -------
        xarray.DataArray
            The geospatial 2D or 3D grid.

        Examples
        --------
        Convert a raster to geospatial and mask it using a Shapely vector geometry:
        sbas.as_geo(grid).rio.clip([geometry])

        Notes
        -----
        This method adds geospatial attributes (CRS and spatial dimensions) to the input grid,
        allowing raster operations using the RioXarray library. If the input grid is already
        in geographic coordinates, the CRS is set to EPSG 4326 with spatial dimensions 'lat' and 'lon'.
        Otherwise, if the input grid is in radar coordinates, a fake metric coordinate system is used
        with EPSG 3857 and spatial dimensions 'y' and 'x'. The method relies on the availability of the
        'rioxarray' module.
        """
        import sys
        assert 'rioxarray' in sys.modules, 'rioxarray module is not found'
        if self.is_geo(da):
            epsg = 4326
            y_dim = 'lat'
            x_dim = 'lon'
        else:
            # fake metrical coordinate system just to perform spatial operations
            epsg = 3857
            y_dim = 'y'
            x_dim = 'x'
        return da.rio.write_crs(epsg).rio.set_spatial_dims(y_dim=y_dim, x_dim=x_dim)

    @staticmethod
    def is_same(grid1, grid2):
        """
        Check if two grids have the same coordinate dimensions.

        Parameters
        ----------
        grid1 : xarray.DataArray
            The first grid to compare.
        grid2 : xarray.DataArray
            The second grid to compare.

        Returns
        -------
        bool
            True if the grids have the same coordinate dimensions, False otherwise.
        """
        dims1 = grid1.dims
        dims2 = grid2.dims
        if 'lat' in dims1 and 'lon' in dims1 and 'lat' in dims2 and 'lon' in dims2:
            return True
        if 'y' in dims1 and 'x' in dims1 and 'y' in dims2 and 'x' in dims2:
            return True
        return False

    def snaphu_config(self, defomax=0, **kwargs):
        """
        Generate a Snaphu configuration file.

        Parameters
        ----------
        defomax : int, optional
            Maximum deformation value. Default is 0.
        **kwargs : dict, optional
            Additional parameters to include in the configuration file.

        Returns
        -------
        str
            The Snaphu configuration file content.

        Notes
        -----
        This method uses the `snaphu_config` method of the PRM object.

        Examples
        --------
        Generate a Snaphu configuration file with defomax=10:
        snaphu_config(defomax=10)

        Generate a Snaphu configuration file with defomax=5 and additional parameters:
        snaphu_config(defomax=5, param1=10, param2=20)
        """
        return self.PRM().snaphu_config(defomax, **kwargs)
 
    @staticmethod
    def cropna(das):
        """
        Crop the valid extent of a raster by removing rows and columns containing only NODATA values.

        Parameters
        ----------
        das : xarray.DataArray
            The input 2D or 3D grid to be cropped.

        Returns
        -------
        xarray.DataArray
            The cropped 2D or 3D grid.

        Examples
        --------
        Crop the valid extent of a raster:
        sbas.cropna(grid)

        Notes
        -----
        This method crops the input grid by removing rows and columns that contain only NODATA values.
        It operates on 2D or 3D grids, where the NODATA values are represented as NaN values.
        The resulting grid has a reduced size, containing only the valid extent of the input grid.
        If the input grid is 3D, the cropping is performed along the dimensions other than 'pair' or 'date'.
        """
        # crop NaNs
        dims = [dim for dim in das.dims if dim != 'pair' and dim != 'date']
        dim0 = [dim for dim in das.dims if dim in ['pair', 'date']]
        assert len(dims) == 2, 'ERROR: the input should be 3D array with "pair" or "date" coordinate'
        da = das.min(dim0)
        indexer = {}
        for dim in dims:
            da = da.dropna(dim=dim, how='all')
            dim_min, dim_max = da[dim].min().item(), da[dim].max().item()
            indexer[dim] = slice(dim_min, dim_max)
        #print ('indexer', indexer)
        return das.loc[indexer]
    
    # replacement for GMTSAR gaussians
    # gauss5x5 = np.genfromtxt('/usr/local/GMTSAR/share/gmtsar/filters/gauss5x5',skip_header=True)
    # gaussian_kernel(5,1) ~= gauss5x5
    @staticmethod
    def gaussian_kernel(size=(5,5), std=(1,1)):
        """
        Generate a 2D Gaussian kernel matrix.

        Parameters
        ----------
        size : tuple, optional
            The size of the kernel matrix in (rows, columns). Default is (5, 5).
        std : tuple, optional
            The standard deviation of the Gaussian distribution in (row_std, column_std). Default is (1, 1).

        Returns
        -------
        numpy.ndarray
            The 2D Gaussian kernel matrix.

        Examples
        --------
        Generate a 5x5 Gaussian kernel with standard deviation of 1 in both dimensions:
        gaussian_kernel(size=(5, 5), std=(1, 1))

        Generate a 3x3 Gaussian kernel with standard deviation of 0.5 in row dimension and 1 in column dimension:
        gaussian_kernel(size=(3, 3), std=(0.5, 1))
        """
        import numpy as np
        from scipy import signal
        matrix1 = signal.gaussian(size[0], std=std[0]).reshape(size[0], 1)
        matrix2 = signal.gaussian(size[1], std=std[1]).reshape(size[1], 1)
        matrix2d = np.outer(matrix1, matrix2)
        return matrix2d

#    @staticmethod
#    def nanconvolve2d(data, kernel, threshold=1/3.):
#        """
#        Convolution using generic kernel on a 2D array with NaNs
#        """
#        import numpy as np
#        import scipy.signal
#        # np.complex128 includes np.float64 real and imagine part
#        vals = data.astype(np.complex128)
#        #vals[np.isnan(data)] = 0
#        #vals[~np.isnan(data)] += 1j
#        nanmask = np.isnan(data)
#        vals[nanmask] = 0
#        vals[~nanmask] += 1j
#
#        conv = scipy.signal.convolve2d(vals, 
#                                      kernel.astype(np.complex128), 
#                                      mode='same', boundary='symm'
#                                     )
#        # suppress incorrect division warning
#        np.seterr(invalid='ignore')
#        return np.where(conv.imag >= threshold*np.sum(kernel), np.divide(conv.real, conv.imag), np.nan)

    @staticmethod
    def nanconvolve2d_gaussian(dataarray, sigmas, truncate):
        """
        Apply lazy convolution using Gaussian kernel on a 2D array with NaN values.

        Parameters
        ----------
        dataarray : xarray.DataArray
            The input 2D array with NaN values.
        sigmas : tuple
            The standard deviations of the Gaussian kernel in (row_sigma, column_sigma).
        truncate : float
            The truncation factor for the Gaussian kernel.

        Returns
        -------
        xarray.DataArray
            The convolved data array with NaN values.

        Examples
        --------
        Apply Gaussian convolution with sigmas (2, 2) and truncate factor 3 to a 2D data array:
        nanconvolve2d_gaussian(dataarray, sigmas=(2, 2), truncate=3)
        """
        import xarray as xr
        import dask
        from dask_image.ndfilters import gaussian_filter as dask_gaussian_filter

        da = dask.array.where(dask.array.isnan(dataarray), 0, 1j + dataarray)
        conv = dask_gaussian_filter(da, sigmas, mode='reflect', truncate=truncate)
        da_conv = xr.DataArray(conv.real/conv.imag, coords=dataarray.coords, name=dataarray.name)
        return da_conv

    def nearest_grid(self, in_grid, search_radius_pixels=None):
        """
        Perform nearest neighbor interpolation on a 2D grid.

        Parameters
        ----------
        in_grid : xarray.DataArray
            The input 2D grid to be interpolated.
        search_radius_pixels : int, optional
            The interpolation distance in pixels. If not provided, the default is set to the chunksize of the SBAS object.

        Returns
        -------
        xarray.DataArray
            The interpolated 2D grid.

        Examples
        --------
        Fill gaps in the specified grid using nearest neighbor interpolation:
        sbas.nearest_grid(grid)

        Notes
        -----
        This method performs nearest neighbor interpolation on a 2D grid. It replaces the NaN values in the input grid with
        the nearest non-NaN values. The interpolation is performed within a specified search radius in pixels.
        If a search radius is not provided, the default search radius is set to the chunksize of the SBAS object.
        """
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np

        assert in_grid.chunks is not None, 'nearest_grid() input grid chunks are not defined'

        if search_radius_pixels is None:
            search_radius_pixels = self.chunksize
        elif search_radius_pixels <= 0:
            print (f'NOTE: interpolation ignored for search_radius_pixels={search_radius_pixels}')
            return in_grid
        else:
            assert search_radius_pixels <= self.chunksize, \
                f'ERROR: apply nearest_grid_pixels() multiple times to fill gaps more than {self.chunksize} pixels chunk size'

        def func(grid, y, x, distance, scaley, scalex):

            grid1d = grid.reshape(-1).copy()
            nanmask0 = np.isnan(grid1d)
            # all the pixels already defined
            if np.all(~nanmask0):
                return grid

            # crop full grid subset to search for missed values neighbors
            ymin = y.min()-scaley*distance-1
            ymax = y.max()+scaley*distance+1
            xmin = x.min()-scalex*distance-1
            xmax = x.max()+scalex*distance+1
            if self.is_ra(in_grid):
                data = in_grid.sel(y=slice(ymin, ymax), x=slice(xmin, xmax))
                ys, xs = data.y, data.x
            else:
                data = in_grid.sel(lat=slice(ymin, ymax), lon=slice(xmin, xmax))
                ys, xs = data.lat, data.lon
            # compute dask arrays to prevent ineffective index lookup
            ys, xs = [vals.values.reshape(-1) for vals in xr.broadcast(ys, xs)]
            data1d = data.values.reshape(-1)
            nanmask = np.isnan(data1d)
            # all the subset pixels are empty, the search is useless
            if np.all(nanmask):
                return grid

            # build index tree for all the valid subset values
            source_yxs = np.stack([ys[~nanmask]/scaley, xs[~nanmask]/scalex], axis=1)
            tree = cKDTree(source_yxs, compact_nodes=False, balanced_tree=False)

            # query the index tree for all missed values neighbors
            target_yxs = np.stack([(y/scaley).reshape(-1)[nanmask0], (x/scalex).reshape(-1)[nanmask0]], axis=1)
            #assert 0, target_yxs
            d, inds = tree.query(target_yxs, k = 1, distance_upper_bound=distance, workers=1)
            # fill missed values using neighbors when these ones are found
            inds = np.where(np.isinf(d), 0, inds)
            grid1d[nanmask0] = np.where(np.isinf(d), np.nan, data1d[~nanmask][inds])
            return grid1d.reshape(grid.shape)

        coords = ['y', 'x'] if self.is_ra(in_grid) else ['lat', 'lon']
        scale = [in_grid[coord].diff(coord).item(0) for coord in coords]
        yy = xr.DataArray(in_grid[coords[0]]).chunk(-1)
        xx = xr.DataArray(in_grid[coords[1]]).chunk(-1)
        ys, xs = xr.broadcast(yy,xx)

        # xarray wrapper
        grid = xr.apply_ufunc(
            func,
            in_grid,
            ys.chunk(in_grid.chunks),
            xs.chunk(in_grid.chunks),
            dask='parallelized',
            vectorize=False,
            output_dtypes=[np.float32],
            dask_gufunc_kwargs={'distance': search_radius_pixels, 'scaley': scale[0], 'scalex': scale[1]},
        )
        assert grid.chunks is not None, 'nearest_grid() output grid chunks are not defined'
        return grid

    def pixel_size(self, grid=(1, 4), average=True):
        """
        Compute the ground pixel size in meters for the default processing grid or the defined one.

        Parameters
        ----------
        grid : tuple or xarray.DataArray, optional
            A pair of x, y grid decimation coefficients or a 2D or 3D Xarray DataArray representing the decimation grid.
            The default is (1, 4) for the default processing grid.
        average : bool, optional
            Flag indicating whether to calculate the average ground pixel size per subswath resolution.
            If True, the average pixel size across all subswaths is returned. If False, a list of pixel sizes per subswath
            is returned. Default is True.

        Returns
        -------
        tuple or list of tuples
            The ground pixel size(s) in meters. If average is True, a tuple of average pixel sizes is returned.
            If average is False, a list of tuples containing the pixel sizes per subswath is returned.

        Examples
        --------
        Get the default average ground pixel size:
        sbas.pixel_size()
        >>> (14.0, 15.7)

        Get the default ground pixel size per subswath:
        sbas.pixel_size(average=False)
        >>> [(14.0, 16.7), (14.0, 14.7)]

        Get the ground pixel size for an unwrapped phase grid with a decimation of {'y': 2, 'x': 2}:
        sbas.pixel_size(unwraps)
        >>> (27.9, 29.5)

        Notes
        -----
        This method computes the ground pixel size in meters for the default processing grid or a user-defined grid.
        The pixel size is calculated based on the azimuth and range pixel sizes obtained from the processing parameters.
        If a grid is provided, the pixel sizes are adjusted according to the grid decimation coefficients.
        By default, the method returns the average pixel size across all subswaths. Setting `average` to False will return
        a list of pixel sizes per subswath. The pixel sizes are rounded to one decimal place.
        """
        import xarray as xr
        import numpy as np

        outs = []
        for subswath in self.get_subswaths():
            # pixel size in meters
            azi_px_size, rng_px_size = self.PRM(subswath).pixel_size()
            # raster pixels decimation
            if isinstance(grid, xr.DataArray):
                dy = grid.y.diff('y')[0].item()
                dx = grid.x.diff('x')[0].item()
            else:
                dy, dx = grid
            outs.append((np.round(azi_px_size*dy,1), np.round(rng_px_size*dx,1)))
        if average:
            pxs = np.asarray(outs)
            return (np.round(pxs[:,0].mean(), 1), np.round(pxs[:,1].mean(), 1))
        else:
            return outs[0] if len(outs) == 1 else outs

    #decimator = lambda da: da.coarsen({'y': 2, 'x': 2}, boundary='trim').mean()
    def pixel_decimator(self, resolution_meters=60, grid=(1, 4), func='mean', debug=False):
        """
        Return function for pixel decimation to the specified output resolution.

        Parameters
        ----------
        resolution_meters : int, optional
            DEM grid resolution in meters. The same grid is used for geocoded results output.
        grid : tuple, optional
            Grid size for pixel decimation in the format (vertical, horizontal).
        debug : bool, optional
            Boolean flag to print debug information.

        Returns
        -------
        callable
            Post-processing function for SBAS.ints() and SBAS.intf_parallel().

        Examples
        --------
        Decimate computed interferograms to default DEM resolution 60 meters:
        decimator = sbas.pixel_decimator()
        sbas.intf_parallel(pairs, func=decimator)
        """
        import numpy as np
        import dask

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

        dy, dx = self.pixel_size(grid)
        yscale, xscale = int(np.round(resolution_meters/dy)), int(np.round(resolution_meters/dx/xscale0))
        if debug:
            print (f'DEBUG: average per subswaths ground pixel size in meters: y={dy}, x={dx}')
        if yscale <= 1 and xscale <= 1 and xscale0==1:
            # decimation impossible
            if debug:
                print (f'DEBUG: decimator = lambda da: da')
            return lambda da: da
        if debug:
            print (f"DEBUG: decimator = lambda da: da.coarsen({{'y': {yscale}, 'x': {xscale0*xscale}}}, boundary='trim').{func}()")

        # decimate function
        def decimator(da):
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

        # return callback function
        return lambda da: decimator(da)

#     # anti-aliasing filter and downscaling, double filter (potentially faster)
#     def antialiasing_downscale(self, da, wavelength, truncate=3, coarsen=(1,4)):
#         import xarray as xr
#         import numpy as np
#         import dask
#         from dask_image.ndfilters import gaussian_filter as dask_gaussian_filter
#         # coarse grid to square cells
#         if coarsen is not None:
#             conv = dask_gaussian_filter(da.data, coarsen, mode='reflect', truncate=2)
#             da_square = xr.DataArray(conv, coords=da.coords, name=da.name).coarsen({'y': coarsen[0], 'x': coarsen[1]}, boundary='trim').mean()
#         # antialiasing (multi-looking) filter
#         if wavelength is None:
#             return da_square
#         dy, dx = self.pixel_size()
#         print ('DEBUG dy, dx', dy, dx)
#         sigmas = int(np.round(wavelength/dy/coarsen[0])), int(np.round(wavelength/dx/coarsen[1]))
#         print ('DEBUG sigmas', sigmas)
#         conv = dask_gaussian_filter(da_square.data, sigmas, mode='reflect', truncate=2)
#         return xr.DataArray(conv, coords=da_square.coords, name=da.name)
    # anti-aliasing filter and downscaling, single filter (potentially more accurate)x
    def antialiasing_downscale(self, da, weight=None, wavelength=None, coarsen=(1,4), debug=False):
        import xarray as xr
        import numpy as np
        import dask
        from dask_image.ndfilters import gaussian_filter as dask_gaussian_filter
        # GMTSAR constant 5.3 defines half-gain at filter_wavelength
        # https://github.com/gmtsar/gmtsar/issues/411
        cutoff = 5.3
        
        # allow this case to save the original grid resolution
        if wavelength is None and (coarsen is None or coarsen==(1,1)):
            return da

        # antialiasing (multi-looking) filter
        if wavelength is None:
            sigmas = np.round([coarsen[0]/cutoff, coarsen[1]/cutoff], 2)
        else:
            dy, dx = self.pixel_size()
            #print ('DEBUG dy, dx', dy, dx)
            #sigmas = int(np.round(wavelength/dy/coarsen[0])), int(np.round(wavelength/dx))
            sigmas = np.round([wavelength/cutoff/dy, wavelength/cutoff/dx], 2)
        if debug:
            print ('DEBUG: antialiasing_downscale sigmas', sigmas, 'for specified wavelength', wavelength)

        # weighted and not weighted convolution
        if weight is None:
            if debug:
                print ('DEBUG: antialiasing_downscale not weighted filtering')
            conv = dask_gaussian_filter(da.data, sigmas, mode='reflect', truncate=2)
        else:
            assert da.shape == weight.shape, f'Different shapes for raster {da.shape} and weight {weight.shape}'
            assert da.dims == weight.dims, f'Different dimension names for raster {da.dims} and weight {weight.dims}'
            if debug:
                print ('DEBUG: antialiasing_downscale weighted filtering')
            conv = dask_gaussian_filter(((1j + da)*weight).data, sigmas, mode='reflect', truncate=2)
            conv = conv.real/conv.imag

        da_conv = xr.DataArray(conv, coords=da.coords, name=da.name)
        # calculate the initial dataarray chunk sizes per dimensions to restore them
        # it works faster when we prevent small chunks usage
        chunksizes = (np.max(da.chunksizes['y']), np.max(da.chunksizes['x']))
        if coarsen is not None:
            # coarse grid to square cells
            return da_conv.coarsen({'y': coarsen[0], 'x': coarsen[1]}, boundary='trim').mean().chunk(chunksizes)
        return da_conv.chunk(chunksizes)
