# ----------------------------------------------------------------------------
# InSAR.dev
# 
# This file is part of the InSAR.dev project: https://InSAR.dev
# 
# Copyright (c) 2025, Alexey Pechnikov
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
    netcdf_engine : str
        The engine used for NetCDF file operations. Default is 'h5netcdf'.
    netcdf_complevel : int
        The compression level for data compression. Default is 3.
    noindex : np.uint32
        The NODATA index value for transform matrices.
    
    Notes
    ----------
    That's possible to define a special NetCDF backend for Docker environment or other cases:

    if os.path.exists('/.dockerenv') and not 'google.colab' in sys.modules:
        # use different NetCDF backend in Docker containers
        from pygmtsar import datagrid
        datagrid.netcdf_engine = 'netcdf4'
    
    """
    import numpy as np

    # Minimum valid Sentinel-1 radar amplitude from GMTSAR code
    #amplitude_threshold = 5.e-21
    # NetCDF options, see https://docs.xarray.dev/en/stable/user-guide/io.html#zarr-compressors-and-filters
    chunksize = 2048
    chunksize1d = 16384
    #netcdf_engine = 'h5netcdf'
    #netcdf_engine = 'netcdf4'
    netcdf_engine_read = 'h5netcdf'
    netcdf_engine_write = 'netcdf4'
    netcdf_format = 'NETCDF4'
    netcdf_chunksize = 512
    netcdf_chunksize1d = 65536
    netcdf_compression_algorithm = 'zlib'
    netcdf_complevel = 3
    netcdf_shuffle = True
    netcdf_queue = 16

    # processing directory
    basedir = '.'

    # define lost class variables due to joblib via arguments
    def _compression(self, shape=None, chunksize=None):
        """
        Return the compression options for a data grid.

        Parameters
        ----------
        shape : tuple, list, np.ndarray, optional
            The shape of the data grid. Required if chunksize is less than grid dimension sizes. Default is None.
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

        if chunksize is None and len(shape) == 1:
            # (stacked) single-dimensional grid 
            chunksize = self.netcdf_chunksize1d
        elif chunksize is None:
            # common 2+D grid
            chunksize = self.netcdf_chunksize

        assert chunksize is not None, 'compression() chunksize is None'
        if isinstance(chunksize, (tuple, list, np.ndarray)):
            # use as is, it can be 2D or 3D grid (even 1D while it is not used for now)
            if shape is not None:
                assert len(shape) == len(chunksize), f'ERROR: defined shape and chunksize dimensions are not equal: {len(shape)} != {len(chunksize)}'
                chunksizes = tuple([chunksize[dim] if chunksize[dim]<shape[dim] else shape[dim] for dim in range(len(shape))])
            else:
                chunksizes = chunksize
        else:
            if shape is not None:
                # 2D or 3D grid
                chunksizes = []
                for idim in range(len(shape)):
                    chunksizes.append(chunksize if chunksize<shape[idim] else shape[idim])
                # set first dimension chunksize to 1 for 3D array
                if len(chunksizes) == 3:
                    chunksizes[0] = 1
                chunksizes = tuple(chunksizes)
            else:
                chunksizes=(chunksize, chunksize)
        opts = dict(chunksizes=chunksizes)
        if self.netcdf_compression_algorithm is not None and self.netcdf_complevel >= 0:
            opts[self.netcdf_compression_algorithm] = True
            opts['complevel'] = self.netcdf_complevel
            opts['shuffle'] = self.netcdf_shuffle
        return opts


    @staticmethod
    def get_spacing(data):
        dy = data.y.diff('y').item(0)
        dx = data.x.diff('x').item(0)
        return (dy, dx)
# 
#     @staticmethod
#     def spatial_ref(da, target=None):
#         """
#         Add geospatial attributes (CRS and spatial dimensions) to allow raster operations using RioXarray.
# 
#         Parameters
#         ----------
#         da : xarray.DataArray or xarray.Dataset
#             The input 2D or 3D grid to be converted to geospatial.
#         target : int, xarray.DataArray, or xarray.Dataset, optional
#             The target EPSG code or an xarray object from which to derive the CRS.
# 
#         Returns
#         -------
#         xarray.DataArray or xarray.Dataset
#             The geospatial 2D or 3D grid with spatial attributes.
# 
#         Examples
#         --------
#         Convert a raster to geospatial and mask it using a Shapely vector geometry:
#         Stack.spatial_ref(grid).rio.clip([geometry])
#         """
#         import xarray as xr
#         import rioxarray
#         import sys
#         #assert 'rioxarray' in sys.modules, 'rioxarray module is not found'
# 
#         if target is None:
#             return da
# 
#         # extract EPSG from target xarray object or use provided EPSG
#         if isinstance(target, (xr.DataArray, xr.Dataset)):
#             if target.rio.crs is None:
#                 raise ValueError("ERROR: Target xarray object has no CRS defined.")
#             epsg = target.rio.crs.to_epsg()
#             if epsg is None:
#                 raise ValueError("ERROR: Target CRS is not an EPSG code, consider using PROJ string.")
#         else:
#             epsg = target
#         #print ('epsg', epsg)
# 
#         if epsg == 4326:
#             # EPSG:4326 (WGS84, lat/lon)
#             da_spatial = (
#                 da.rio.write_crs(4326)
#                   .rio.set_spatial_dims(y_dim='lat', x_dim='lon')
#                   .rio.write_grid_mapping()
#                   .assign_coords(lat=da.lat.assign_attrs(axis='Y', 
#                                                        standard_name='latitude',
#                                                        long_name='latitude'),
#                                lon=da.lon.assign_attrs(axis='X', 
#                                                        standard_name='longitude',
#                                                        long_name="longitude"))
#                   .assign_attrs({'Conventions': 'CF-1.8'}))
#             if isinstance(da, xr.Dataset):
#                 return da_spatial.assign({
#                     var: da[var].assign_attrs(grid_mapping='spatial_ref', coordinates='lon lat')  
#                         for var in da.data_vars if set(da[var].dims) == {'lat', 'lon'}
#                 })
#             return da_spatial.assign_attrs(grid_mapping='spatial_ref', coordinates='lon lat')
# 
#         # projected coordinates
#         da_spatial = (
#             da.rio.write_crs(epsg)
#               .rio.set_spatial_dims(y_dim='y', x_dim='x')
#               .rio.write_grid_mapping()
#               .assign_coords(y=da.y.assign_attrs(axis='Y', 
#                                                standard_name='projection_y_coordinate',
#                                                long_name='northing'),
#                            x=da.x.assign_attrs(axis='X', 
#                                                standard_name='projection_x_coordinate',
#                                                long_name="easting"))
#               .assign_attrs({'Conventions': 'CF-1.8'}))
# 
#         if isinstance(da, xr.Dataset):
#             return da_spatial.assign({
#                 var: da[var].assign_attrs(grid_mapping='spatial_ref', coordinates='x y')  
#                     for var in da.data_vars if set(da[var].dims) == {'y', 'x'}
#             })
#         return da_spatial.assign_attrs(grid_mapping='spatial_ref', coordinates='x y')

    @staticmethod
    def spatial_ref(da, target=None):
        """
        Add geospatial attributes (CRS and spatial dimensions) to allow raster operations using RioXarray.

        Parameters
        ----------
        da : xarray.DataArray or xarray.Dataset
            The input 2D or 3D grid to be converted to geospatial.
        target : int, xarray.DataArray, or xarray.Dataset, optional
            The target EPSG code or an xarray object from which to derive the CRS.

        Returns
        -------
        xarray.DataArray or xarray.Dataset
            The geospatial 2D or 3D grid with spatial attributes.

        Examples
        --------
        Convert a raster to geospatial and mask it using a Shapely vector geometry:
        Stack.spatial_ref(grid).rio.clip([geometry])
        """
        import xarray as xr
        import rioxarray
        import sys
        #assert 'rioxarray' in sys.modules, 'rioxarray module is not found'

        if target is None:
            return da

        # extract EPSG from target xarray object or use provided EPSG
        if isinstance(target, (xr.DataArray, xr.Dataset)):
            if target.rio.crs is None:
                raise ValueError("ERROR: Target xarray object has no CRS defined.")
            epsg = target.rio.crs.to_epsg()
            if epsg is None:
                raise ValueError("ERROR: Target CRS is not an EPSG code, consider using PROJ string.")
        else:
            epsg = target
        #print ('spatial_ref epsg', epsg)

        if epsg == 4326:
            # EPSG:4326 (WGS84, lat/lon)
            da_spatial = (
                da.rio.write_crs(4326)
                  .rio.set_spatial_dims(y_dim='lat', x_dim='lon')
                  #.rio.write_grid_mapping()
                  .assign_coords(lat=da.lat.assign_attrs(axis='Y', 
                                                       standard_name='latitude',
                                                       long_name='latitude'),
                               lon=da.lon.assign_attrs(axis='X', 
                                                       standard_name='longitude',
                                                       long_name="longitude"))
                  .assign_attrs({'Conventions': 'CF-1.8'}))
            if isinstance(da_spatial, xr.Dataset):
                da_spatial = da_spatial.assign({
                    var: da_spatial[var].assign_attrs(grid_mapping='spatial_ref', coordinates='lon lat')  
                        for var in da_spatial.data_vars if 'lat' in da_spatial[var].dims and 'lon' in da_spatial[var].dims
                })
            else:
                da_spatial = da_spatial.assign_attrs(grid_mapping='spatial_ref', coordinates='lon lat')
        else:
            # projected coordinates
            da_spatial = (
                da.rio.write_crs(epsg)
                  .rio.set_spatial_dims(y_dim='y', x_dim='x')
                  #.rio.write_grid_mapping()
                  .assign_coords(y=da.y.assign_attrs(axis='Y', 
                                                   standard_name='projection_y_coordinate',
                                                   long_name='northing'),
                               x=da.x.assign_attrs(axis='X', 
                                                   standard_name='projection_x_coordinate',
                                                   long_name="easting"))
                  .assign_attrs({'Conventions': 'CF-1.8'}))
    
            if isinstance(da_spatial, xr.Dataset):
                da_spatial = da_spatial.assign({
                    var: da_spatial[var].assign_attrs(grid_mapping='spatial_ref', coordinates='x y')  
                        for var in da_spatial.data_vars if 'y' in da_spatial[var].dims and 'x' in da_spatial[var].dims
                })
            else:
                da_spatial = da_spatial.assign_attrs(grid_mapping='spatial_ref', coordinates='x y')

        return da_spatial.assign_attrs(coordinates='spatial_ref')

#     # da.dropna(dim=dim, how='all') is not fast at all
#     @staticmethod
#     def cropna(das, index=-1):
#         """
#         Crop the valid extent of a raster by removing rows and columns containing only NODATA values.
# 
#         Parameters
#         ----------
#         das : xarray.DataArray
#             The input 2D or 3D grid to be cropped.
# 
#         Returns
#         -------
#         xarray.DataArray
#             The cropped 2D or 3D grid.
# 
#         Examples
#         --------
#         Crop the valid extent of a raster:
#         stack.cropna(grid)
# 
#         Notes
#         -----
#         This method crops the input grid by removing rows and columns that contain only NODATA values.
#         It operates on 2D or 3D grids, where the NODATA values are represented as NaN values.
#         The resulting grid has a reduced size, containing only the valid extent of the input grid.
#         If the input grid is 3D, the cropping is performed along the dimensions other than 'pair' or 'date'.
#         """
#         # crop NaNs
#         dims = [dim for dim in das.dims if dim != 'pair' and dim != 'date']
#         dim0 = [dim for dim in das.dims if dim in ['pair', 'date']]
#         #print ('dims', dims, 'dim0', dim0)
#         assert len(dims) == 2, 'ERROR: the input should be 3D array with "pair" or "date" coordinate'
#         # slow check using all the grids in the stack
#         #da = das.min(dim0)
#         # fast check using the only "index" grid in the stack
#         da = das.isel({dim0[index]: index}) if dim0 != [] else das
#         indexer = {}
#         for dim in dims:
#             da = da.dropna(dim=dim, how='all')
#             dim_min, dim_max = da[dim].min().item(), da[dim].max().item()
#             indexer[dim] = slice(dim_min, dim_max)
#         #print ('indexer', indexer)
#         return das.loc[indexer]

#     # replacement for GMTSAR gaussians
#     # gauss5x5 = np.genfromtxt('/usr/local/GMTSAR/share/gmtsar/filters/gauss5x5',skip_header=True)
#     # gaussian_kernel(5,1) ~= gauss5x5
#     @staticmethod
#     def gaussian_kernel(size=(5,5), std=(1,1)):
#         """
#         Generate a 2D Gaussian kernel matrix.
# 
#         Parameters
#         ----------
#         size : tuple, optional
#             The size of the kernel matrix in (rows, columns). Default is (5, 5).
#         std : tuple, optional
#             The standard deviation of the Gaussian distribution in (row_std, column_std). Default is (1, 1).
# 
#         Returns
#         -------
#         numpy.ndarray
#             The 2D Gaussian kernel matrix.
# 
#         Examples
#         --------
#         Generate a 5x5 Gaussian kernel with standard deviation of 1 in both dimensions:
#         gaussian_kernel(size=(5, 5), std=(1, 1))
# 
#         Generate a 3x3 Gaussian kernel with standard deviation of 0.5 in row dimension and 1 in column dimension:
#         gaussian_kernel(size=(3, 3), std=(0.5, 1))
#         """
#         import numpy as np
#         from scipy import signal
#         matrix1 = signal.gaussian(size[0], std=std[0]).reshape(size[0], 1)
#         matrix2 = signal.gaussian(size[1], std=std[1]).reshape(size[1], 1)
#         matrix2d = np.outer(matrix1, matrix2)
#         return matrix2d

    @staticmethod
    def get_bounds(geometry):
        import geopandas as gpd
        import xarray as xr
    
        if isinstance(geometry, (xr.DataArray, xr.Dataset)) and ('lat' in geometry.dims and 'lon' in geometry.dims):
            lon_start = geometry.lon.min().item()
            lat_start = geometry.lat.min().item()
            lon_end   = geometry.lon.max().item()
            lat_end   = geometry.lat.max().item()
            bounds = lon_start, lat_start, lon_end, lat_end
        elif isinstance(geometry, (xr.DataArray, xr.Dataset)):
            x_start = geometry.x.min().item()
            y_start = geometry.y.min().item()
            x_end   = geometry.x.max().item()
            y_end   = geometry.y.max().item()
            bounds = x_start, y_start, x_end, y_end
        elif isinstance(geometry, gpd.GeoDataFrame):
            bounds = geometry.dissolve().envelope.item().bounds
        elif isinstance(geometry, gpd.GeoSeries):
            bounds = geometry.union_all().envelope.bounds
        elif isinstance(geometry, tuple):
            # geometry is already bounds
            bounds = geometry
        else:
            bounds = geometry.bounds
        #print ('bounds', bounds)
        #lon_start, lat_start, lon_end, lat_end
        return bounds

    @staticmethod
    def calculate_coarsen_start(da, name, spacing, grid_factor=1):
        """
        Calculate start coordinate to align coarsened grids.
        """
        import numpy as np
        for i in range(spacing):
            values = da[name].isel({name: slice(i, None)}).coarsen({name: spacing}, boundary='trim').mean().values
            delta = np.floor(values[0] % (spacing*grid_factor))
            #print ('i', i, 'delta', delta, 'values', values[:5])
            if delta == 0:
                #print ('calculate_start', name, i)
                return i
        return None

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

    def nearest_grid(self, in_grid, search_radius_pixels=None):
        """
        Perform nearest neighbor interpolation on a 2D grid.

        Parameters
        ----------
        in_grid : xarray.DataArray
            The input 2D grid to be interpolated.
        search_radius_pixels : int, optional
            The interpolation distance in pixels. If not provided, the default is set to the chunksize of the Stack object.

        Returns
        -------
        xarray.DataArray
            The interpolated 2D grid.

        Examples
        --------
        Fill gaps in the specified grid using nearest neighbor interpolation:
        stack.nearest_grid(grid)

        Notes
        -----
        This method performs nearest neighbor interpolation on a 2D grid. It replaces the NaN values in the input grid with
        the nearest non-NaN values. The interpolation is performed within a specified search radius in pixels.
        If a search radius is not provided, the default search radius is set to the chunksize of the Stack object.
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
            data = in_grid.sel(y=slice(ymin, ymax), x=slice(xmin, xmax))
            ys, xs = data.y, data.x
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

        coords = ['y', 'x']
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
