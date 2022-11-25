#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar

class datagrid:
    import numpy as np

    # NetCDF options, see https://docs.xarray.dev/en/stable/user-guide/io.html#zarr-compressors-and-filters
    chunksize = 512
    compression = dict(zlib=True, complevel=3, chunksizes=(chunksize, chunksize))
    engine = 'h5netcdf'
    # NODATA index value for transform matrices
    # TODO: use the same datatype for the matrices to allow 64 bit datatype
    noindex = np.uint32(-1)

    @staticmethod
    def is_ra(grid):
        dims = grid.dims
        if 'y' in dims and 'x' in dims:
            return True
        return False

    @staticmethod
    def is_geo(grid):
        dims = grid.dims
        if 'lat' in dims and 'lon' in dims:
            return True
        return False

    def as_geo(self, da):
        """
        Add spatial attributes to allow use rioxarray functions like to .rio.clip([geometry])
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
        dims1 = grid1.dims
        dims2 = grid2.dims
        if 'lat' in dims1 and 'lon' in dims1 and 'lat' in dims2 and 'lon' in dims2:
            return True
        if 'y' in dims1 and 'x' in dims1 and 'y' in dims2 and 'x' in dims2:
            return True
        return False

    def snaphu_config(self, defomax=0, **kwargs):
        return self.PRM().snaphu_config(defomax, **kwargs)
 
    @staticmethod
    def cropna(das):
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
        """Make 2D Gaussian kernel matrix"""
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
        Lazy convolution using Gaussian kernel on a 2D array with NaNs
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
        Nearest Neighbour interpolation
        """
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np

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
        return grid

    def pixel_size(self, grid=(1, 4), average=True):
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
    def pixel_decimator(self, resolution_meters=60, grid=(1, 4), debug=False):
        import numpy as np

        dy, dx = self.pixel_size(grid)
        yy, xx = int(np.round(resolution_meters/dy)), int(np.round(resolution_meters/dx))
        if debug:
            print (f'DEBUG: average per subswaths ground pixel size in meters: y={dy}, x={dx}')
        if yy <= 1 and xx <= 1:
            if debug:
                print (f"DEBUG: decimator = lambda da: da")
            return lambda da: da
        if debug:
            print (f"DEBUG: decimator = lambda da: da.coarsen({{'y': {yy}, 'x': {xx}}}, boundary='trim').mean()")
        return lambda da: da.coarsen({'y': yy, 'x': xx}, boundary='trim').mean()
