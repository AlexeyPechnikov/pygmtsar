#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar

class datagrid:
    import numpy as np

    # NetCDF options
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

    @staticmethod
    def nearest_grid(in_grid, search_radius_pixels=300):
        """
        Pixel-based Nearest Neighbour interpolation
        """
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np
        ys, xs = np.meshgrid(range(in_grid.shape[1]), range(in_grid.shape[0]))
        ys = ys.reshape(-1)
        xs = xs.reshape(-1)
        if isinstance(in_grid, xr.DataArray):
            zs = in_grid.values.reshape(-1)
        else:
            zs = in_grid.reshape(-1)
        mask = np.where(~np.isnan(zs))
        # on regular source grid some options should be redefined for better performance
        tree = cKDTree(np.column_stack((ys[mask],xs[mask])), compact_nodes=False, balanced_tree=False)
        # use distance_limit
        d, inds = tree.query(np.column_stack((ys,xs)), k = 1, distance_upper_bound=search_radius_pixels, workers=8)
        # replace not available indexes by zero (see distance_upper_bound)
        fakeinds = np.where(~np.isinf(d), inds, 0)
        # produce the same output array as dataset to be able to add global attributes
        values = np.where(~np.isinf(d), zs[mask][fakeinds], np.nan).reshape(in_grid.shape)
        if isinstance(in_grid, xr.DataArray):
            return xr.DataArray(values, coords=in_grid.coords, name=in_grid.name)
        return values

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
            return (pxs[:,0].mean(), pxs[:,1].mean())
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
