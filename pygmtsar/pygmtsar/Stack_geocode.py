# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_sbas import Stack_sbas
from .tqdm_dask import tqdm_dask

class Stack_geocode(Stack_sbas):

    def geocode(self, coarsen=60.):
        """
        Build topography in radar coordinates from WGS84 DEM using parallel computation.

        Parameters
        ----------
        interactive : bool, optional
            If True, the computation will be performed interactively and the results will be returned as delayed objects.
            If False, the progress will be displayed using tqdm_dask. Default is False.

        Returns
        -------
        handler or list of handlers
            The handler(s) of the delayed computation if 'interactive' is True. Otherwise, None.

        Examples
        --------
        stack.topo()

        Notes
        -----
        This method performs the parallel computation of topography in the radar coordinates using Dask.
        If 'interactive' is True, the delayed computation handlers will be returned.
        Otherwise, the progress will be displayed using tqdm_dask.
        """
        import warnings
        warnings.filterwarnings('ignore')

        self.trans(coarsen=coarsen)
        self.trans_inv(coarsen=coarsen)

#     def stack_ra2ll(self, data, **kwargs):
#         return self.ra2ll(data, **kwargs)
# 
#     def cube_ra2ll(self, data, **kwargs):
#         return self.ra2ll(data, **kwargs)

    # coarsen=1:
    # nearest: coords [array([596.97436523]), array([16976.93164062])]
    # linear:  coords [array([597.10408125]), array([16977.41447869])]
    # coarsen=2:
    # nearest: coords [array([596.4810791]), array([16977.52734375])]
    # linear:  coords [array([597.10471665]), array([16977.3737493])]
    # coarsen=4:
    # nearest: coords [array([596.42352295]), array([16978.65625])]
    # linear:  coords [array([597.1080563]), array([16977.35608873])]
    def geocode_ll2ra(self, data, z_offset=None):
        """
        Inverse geocode input geodataframe with 2D or 3D points. 
        """
        import numpy as np
        import geopandas as gpd
        from shapely import Point
        from scipy.interpolate import RegularGridInterpolator

        dem = self.get_dem()
        prm = self.PRM()

        dy = dem.lat.diff('lat')[0]
        dx = dem.lon.diff('lon')[0]
        points_ll = []
        for geom in data.geometry:
            subset = dem.sel(lat=slice(geom.y-2*dy, geom.y+2*dy),lon=slice(geom.x-2*dx, geom.x+2*dx)).compute(n_process=1)
            #print (subset.shape)
            # perform interpolation
            lats, lons = subset.lat.values, subset.lon.values
            interp = RegularGridInterpolator((lats, lons),
                                             subset.values.astype(np.float64),
                                             method='cubic',
                                             bounds_error=False)
            # interpolate specified point elevation on DEM adding 3D point vertical coordinate when exists
            if z_offset is None:
                ele = interp([geom.y, geom.x])[0] + (geom.z if geom.has_z else 0)
            else:
                ele = interp([geom.y, geom.x])[0] + z_offset
            points_ll.append([geom.x, geom.y, ele])
            del interp
        points_ra = prm.SAT_llt2rat(points_ll)
        points_ra = points_ra[:,:2] if len(points_ll)>1 else [points_ra[:2]]
        #print ('points_ra', points_ra)
        # set fake CRS to differ from WGS84 coordinates
        return gpd.GeoDataFrame(data, geometry=[Point(point_ra) for point_ra in points_ra], crs=3857)

##########################################################################################
# ra2ll
##########################################################################################
    def ra2ll(self, data, autoscale=True):
        """
        Perform geocoding from radar to geographic coordinates.

        Parameters
        ----------
        grid : xarray.DataArray
            Grid(s) representing the interferogram(s) in radar coordinates.
        trans : xarray.DataArray
            Geocoding transform matrix in radar coordinates.

        Returns
        -------
        xarray.DataArray
            The inverse geocoded grid(s) in geographic coordinates.

        Examples
        --------
        Geocode 3D unwrapped phase grid stack:
        unwraps_ll = stack.intf_ra2ll(stack.open_grids(pairs, 'unwrap'))
        # or use "geocode" option for open_grids() instead:
        unwraps_ll = stack.open_grids(pairs, 'unwrap', geocode=True)
        """
        import dask
        import xarray as xr
        import numpy as np

        # helper check
        if not 'y' in data.dims or not 'x' in data.dims:
            print ('NOTE: the input data not in radar coordinates, miss geocoding')
            return data

        # get complete transform table
        trans = self.get_trans()

        @dask.delayed
        def intf_block(lats_block, lons_block, stackval=None):
            from scipy.interpolate import RegularGridInterpolator
        
            # use outer variables
            block_grid = (data.sel({stackvar: stackval}) if stackval is not None else data).compute(n_workers=1)
            trans_block = trans.sel(lat=lats_block, lon=lons_block).compute(n_workers=1)
        
            # check if the data block exists
            if not (trans_block.lat.size>0 and trans_block.lon.size>0):
                del block_grid, trans_block
                return np.nan * np.zeros((lats_block.size, lons_block.size), dtype=np.float32)

            # use trans table subset
            y = trans_block.azi.values.ravel()
            x = trans_block.rng.values.ravel()
            points = np.column_stack([y, x])
            del trans_block
        
            # get interferogram full grid
            ys = data.y.values
            xs = data.x.values

            # this code spends additional time for the checks to exclude warnings
            if np.all(np.isnan(y)):
                del block_grid, ys, xs, points
                return np.nan * np.zeros((lats_block.size, lons_block.size), dtype=np.float32)

            # calculate trans grid subset extent
            ymin, ymax = np.nanmin(y), np.nanmax(y)
            xmin, xmax = np.nanmin(x), np.nanmax(x)
            del y, x
            # and spacing
            dy = ys[1] - ys[0]
            dx = xs[1] - xs[0]
        
            # select required interferogram grid subset
            ys_subset = ys[(ys>ymin-dy)&(ys<ymax+dy)]
            xs_subset = xs[(xs>xmin-dx)&(xs<xmax+dx)]
            del ymin, ymax, xmin, xmax, dy, dx, ys, xs
        
            # for cropped interferogram we can have no valid pixels for the processing
            if ys_subset.size == 0 or xs_subset.size == 0:
                del ys_subset, xs_subset, points, block_grid
                return np.nan * np.zeros((lats_block.size, lons_block.size), dtype=np.float32)

            # Wall time: 1min 25s

            values = block_grid.sel(y=ys_subset, x=xs_subset).compute(n_workers=1).data.astype(np.float64)
            del block_grid
        
            # Wall time: 7min 47s
            # distributed.worker.memory - WARNING - Worker is at 81% memory usage.

            # perform interpolation
            interp = RegularGridInterpolator((ys_subset, xs_subset), values, method='nearest', bounds_error=False)
            grid_ll = interp(points).reshape(lats_block.size, lons_block.size).astype(np.float32)
            del ys_subset, xs_subset, points, values
            return grid_ll

        # analyse grid and transform matrix spacing
        grid_dy = np.diff(data.y)[0]
        grid_dx = np.diff(data.x)[0]
        trans_dy = np.diff(trans.y)[0]
        trans_dx = np.diff(trans.x)[0]

        # define transform spacing in radar coordinates
        step_y = int(np.round(grid_dy / trans_dy))
        step_x = int(np.round(grid_dx / trans_dx))
        #print ('step_y', step_y, 'step_x', step_x)
        assert step_y>=1 and step_x>=1, f'Transforming grid spacing (grid_dy, grid_dx) is smaller \
                                          than transform matrix spacing (trans_dy, trans_dx), \
                                          call Stack.geocode() with less coarsing'
        # decimate the full trans grid to the required spacing
        if autoscale and (step_y>1 or step_x>1):
            # define the equally spacing geographic coordinates grid
            trans = trans.sel(lat=trans.lat[::step_y], lon=trans.lon[::step_x])
         # define output geographic coordinates grid
        lats = trans.lat
        lons = trans.lon

        # select required variables only
        trans = trans[['azi', 'rng']]

        # split to equal chunks and rest
        lats_blocks = np.array_split(lats, np.arange(0, lats.size, self.chunksize)[1:])
        lons_blocks = np.array_split(lons, np.arange(0, lons.size, self.chunksize)[1:])

        # find stack dim
        stackvar = data.dims[0] if len(data.dims) == 3 else 'stack'
        #print ('stackvar', stackvar)

        stack = []
        for stackval in data[stackvar].values if len(data.dims) == 3 else [None]:
            # per-block processing
            blocks_total  = []
            for lats_block in lats_blocks:
                blocks = []
                for lons_block in lons_blocks:
                    block = dask.array.from_delayed(intf_block(lats_block, lons_block, stackval),
                                                    shape=(lats_block.size, lons_block.size), dtype=np.float32)
                    blocks.append(block)
                    del block
                blocks_total.append(blocks)
                del blocks
            dask_block = dask.array.block(blocks_total)
            if len(data.dims) == 3:
                coords = {stackvar: [stackval], 'lat': trans.coords['lat'], 'lon': trans.coords['lon']}
                da = xr.DataArray(dask_block[None, :], coords=coords)
            else:
                coords = {'lat': trans.coords['lat'], 'lon': trans.coords['lon']}
                da = xr.DataArray(dask_block, coords=coords)
            stack.append(da)
            del blocks_total, dask_block

        # wrap lazy Dask array to Xarray dataarray
        if len(data.dims) == 2:
            out = stack[0]
        else:
            out = xr.concat(stack, dim=stackvar)
        del stack

        # append source grid coordinates excluding removed y, x ones
        for (k,v) in data.coords.items():
            if k not in ['y','x']:
                out[k] = v
        return out.rename(data.name)

##########################################################################################
# ll2ra
##########################################################################################
    def ll2ra(self, data):
        """
        Perform geocoding from geographic to radar coordinates.

        Parameters
        ----------
        grid : xarray.DataArray
            Grid(s) in geographic coordinates.
        trans : xarray.DataArray
            Inverse geocoding transform matrix.

        Returns
        -------
        xarray.DataArray
            The inverse geocoded grid(s) in radar coordinates.

        Examples
        --------
        Inverse geocode 2D land mask grid:
        landmask_ra = stack.ll2ra(stack.get_landmask())
        """
        import dask
        import xarray as xr
        import numpy as np

        # helper check
        if not 'lat' in data.dims or not 'lon' in data.dims:
            print ('NOTE: the input data not in geograophic coordinates, miss inverse geocoding')
            return data

        # get complete inverse transform table
        trans_inv = self.get_trans_inv()

        @dask.delayed
        def geo_block(ys_block, xs_block, stackval=None):
            from scipy.interpolate import RegularGridInterpolator

            def nangrid():
                return np.nan * np.zeros((ys_block.size, xs_block.size), dtype=np.float32)
        
            # use outer variables
            block_grid = data.sel({stackvar: stackval}) if stackval is not None else data
            trans_inv_block = trans_inv.sel(y=ys_block, x=xs_block).compute(n_workers=1)

            # check if the data block exists
            if not (trans_inv_block.y.size>0 and trans_inv_block.x.size>0):
                del block_grid, trans_inv_block
                return nangrid()

            # use trans table subset
            y = trans_inv_block.lt.values.ravel()
            x = trans_inv_block.ll.values.ravel()
            points = np.column_stack([y, x])
            del trans_inv_block

            # get interferogram full grid
            ys = data.lat.values
            xs = data.lon.values

            # this code spends additional time for the checks to exclude warnings
            if np.all(np.isnan(y)):
                del block_grid, ys, xs, points
                return nangrid()

            # calculate trans grid subset extent
            ymin, ymax = np.nanmin(y), np.nanmax(y)
            xmin, xmax = np.nanmin(x), np.nanmax(x)
            del y, x
            # and spacing
            dy = ys[1] - ys[0]
            dx = xs[1] - xs[0]

            # select required interferogram grid subset
            ys_subset = ys[(ys>ymin-dy)&(ys<ymax+dy)]
            xs_subset = xs[(xs>xmin-dx)&(xs<xmax+dx)]
            del ymin, ymax, xmin, xmax, dy, dx, ys, xs

            # for cropped interferogram we can have no valid pixels for the processing
            if ys_subset.size == 0 or xs_subset.size == 0:
                del ys_subset, xs_subset, points, block_grid
                return nangrid()

            # Wall time: 1min 25s

            values = block_grid.sel(lat=ys_subset, lon=xs_subset).compute(n_workers=1).data.astype(np.float64)
            del block_grid

            # Wall time: 7min 47s
            # distributed.worker.memory - WARNING - Worker is at 81% memory usage.

            # perform interpolation
            interp = RegularGridInterpolator((ys_subset, xs_subset), values, method='nearest', bounds_error=False)
            grid_ra = interp(points).reshape(ys_block.size, xs_block.size).astype(np.float32)
            del ys_subset, xs_subset, points, values
            return grid_ra

         # define output geographic coordinates grid
        lats = trans_inv.y
        lons = trans_inv.x

        # split to equal chunks and rest
        lats_blocks = np.array_split(lats, np.arange(0, lats.size, self.chunksize)[1:])
        lons_blocks = np.array_split(lons, np.arange(0, lons.size, self.chunksize)[1:])

        # find stack dim
        stackvar = data.dims[0] if len(data.dims) == 3 else 'stack'
        #print ('stackvar', stackvar)

        stack = []
        for stackval in data[stackvar].values if len(data.dims) == 3 else [None]:
            # per-block processing
            blocks_total  = []
            for lats_block in lats_blocks:
                blocks = []
                for lons_block in lons_blocks:
                    block = dask.array.from_delayed(geo_block(lats_block, lons_block, stackval),
                                                    shape=(lats_block.size, lons_block.size), dtype=np.float32)
                    blocks.append(block)
                    del block
                blocks_total.append(blocks)
                del blocks
            dask_block = dask.array.block(blocks_total)
            if len(data.dims) == 3:
                coords = {stackvar: [stackval], 'y': trans_inv.coords['y'], 'x': trans_inv.coords['x']}
                da = xr.DataArray(dask_block[None, :], coords=coords)
            else:
                coords = {'y': trans_inv.coords['y'], 'x': trans_inv.coords['x']}
                da = xr.DataArray(dask_block, coords=coords)
            stack.append(da)
            del blocks_total, dask_block

        # wrap lazy Dask array to Xarray dataarray
        if len(data.dims) == 2:
            out = stack[0]
        else:
            out = xr.concat(stack, dim=stackvar)
        del stack

        # append source grid coordinates excluding removed y, x ones
        for (k,v) in data.coords.items():
            if k not in ['lat','lon']:
                out[k] = v
        return out.rename(data.name)

#     def intf_ll2ra_matrix(self, stack, chunksize=None, interactive=False):
#         """
#         Perform parallel computation of the geographic-to-radar coordinate transformation matrix.
# 
#         Parameters
#         ----------
#         pairs : list or None, optional
#             List of interferogram pairs to process. If None, all available pairs will be processed.
#         interactive : bool, optional
#             Flag indicating whether to return the matrix without saving to a file.
#         """
#         import xarray as xr
#         import numpy as np
#         import dask
#         import os
# 
#         if len(stack.dims) == 3:
#             intf = stack[0].astype(bool)
#         else:
#             intf = stack.astype(bool)
# 
#         if chunksize is None:
#             chunksize = self.chunksize
# 
#         @dask.delayed
#         def intf_block_inv(intf_block):
#             from scipy.interpolate import RegularGridInterpolator
# 
#             # define input grid subset
#             trans_inv_dy = np.diff(trans_inv.y)[0]
#             trans_inv_dx = np.diff(trans_inv.x)[0]
#             ymin, ymax = np.min(intf_block.y)-trans_inv_dy, np.max(intf_block.y)+trans_inv_dy
#             xmin, xmax = np.min(intf_block.x)-trans_inv_dx, np.max(intf_block.x)+trans_inv_dx
#             trans_block = trans_inv.sel(y=slice(ymin, ymax), x=slice(xmin,xmax))
#             #print ('trans_block', trans_block)
# 
#             # for cropped interferogram we can have no valid pixels for the processing
#             if trans_block.y.size == 0 or trans_block.x.size == 0:
#                 fake_block = np.nan * np.zeros((intf.y.size, intf.x.size), dtype=np.float32)
#                 return np.asarray([fake_block, fake_block])
# 
#             # define input grid
#             trans_grid = (trans_block.y, trans_block.x)
# 
#             # define output grid
#             ys, xs = np.meshgrid(intf_block.y, intf_block.x)
#             points = np.column_stack([ys.ravel(), xs.ravel()])
# 
#             grids = []
#             for var in ['lt', 'll']:
#                 # np.float64 required for the interpolator
#                 block = trans_block[var].compute(n_workers=1).data.astype(np.float64)
#                 # perform interpolation
#                 interp = RegularGridInterpolator(trans_grid, block, method='linear', bounds_error=False)
#                 grid = interp(points).reshape((intf_block.y.size, intf_block.x.size))
#                 grids.append(grid)
#                 del block, interp, grid
# 
#             del trans_block, trans_grid, ys, xs, points
#             return np.asarray(grids)
# 
#         blocks = chunksize**2 // intf.y.size
#         xs_blocks = np.array_split(intf.x, np.arange(0, intf.x.size, blocks if blocks>=1 else 1)[1:])
#         #print ('xs_blocks', xs_blocks)
# 
#         # get transform table
#         trans_inv = self.get_trans_dat_inv()[['lt', 'll']]
#         #print ('trans_inv', trans_inv)
# 
#         # per-block processing
#         block_lts  = []
#         block_lls  = []
#         for xs_block in xs_blocks:
#             block_lt, block_ll = dask.array.from_delayed(intf_block_inv(intf.sel(x=xs_block)),
#                                             shape=(2, intf.y.size, xs_block.size), dtype=np.float32)
#             block_lts.append(block_lt)
#             block_lls.append(block_ll)
#             del block_lt, block_ll
#         # set the output grid and drop the fake dimension if needed
#         lt = xr.DataArray(dask.array.block(block_lts), coords=intf.coords)
#         ll = xr.DataArray(dask.array.block(block_lls), coords=intf.coords)
#         trans = xr.Dataset({'lt': lt, 'll': ll}).rename({'y': 'a', 'x': 'r'})
# 
#         if interactive:
#             return trans
#         return self.save_grid(trans, 'intf_ll2ra', f'Build ll2ra Transform', chunksize)
# 
#     def get_intf_ll2ra(self, chunksize=None):
#         """
#         Get the geographic-to-radar coordinate transformation matrix.
# 
#         Parameters
#         ----------
#         chunksize : int or dict, optional
#             Chunk size for dask arrays. If not provided, the default chunk size is used.
# 
#         Returns
#         -------
#         xarray.DataArray
#             The radar-to-geographic coordinate transformation matrix.
# 
#         Notes
#         -----
#         This method retrieves the geographic-to-radar coordinate transformation matrix (intf_ll2ra) stored in the
#         NetCDF grid. The matrix is useful for direct geocoding, converting geographic coordinate grids to
#         radar coordinates interferogram grid.
#         """
#         import xarray as xr
# 
#         filename = self.get_filenames(None, 'intf_ll2ra')
#         trans_inv = xr.open_dataset(filename, engine=self.engine, chunks=self.chunksize)
#         if 'a' in trans_inv.dims and 'r' in trans_inv.dims:
#             return trans_inv.rename({'a': 'y', 'r': 'x'})
#         else:
#             return trans_inv
# 
#     def intf_ll2ra(self, grids, chunksize=None):
#         """
#         Perform inverse geocoding from geographic to radar coordinates.
# 
#         Parameters
#         ----------
#         grids : xarray.DataArray
#             Grid(s) representing the interferogram(s) in radar coordinates.
#         chunksize : int or dict, optional
#             Chunk size for dask arrays. If not provided, the default chunk size is used.
# 
#         Returns
#         -------
#         xarray.DataArray
#             The inverse geocoded grid(s) in geographic coordinates.
# 
#         Examples
#         --------
#         Inverse geocode 3D unwrapped phase grids stack:
#         unwraps_ll = stack.open_grids(pairs, 'unwrap', geocode=True)
#         unwraps = stack.intf_ll2ra(unwraps_ll)
#         """
#         import dask
#         import xarray as xr
#         import numpy as np
# 
#         # helper check
#         if not 'lat' in grids.dims or not 'lon' in grids.dims:
#             print ('NOTE: the input grid is not in geographic coordinates, miss inverse geocoding')
#             return grids
# 
#         if chunksize is None:
#             chunksize = self.chunksize
# 
#         @dask.delayed
#         def intf_block_inv(trans_inv, grid_ll):
#             from scipy.interpolate import RegularGridInterpolator
# 
#             # use trans table subset
#             lt = trans_inv.lt.values.ravel()
#             ll = trans_inv.ll.values.ravel()
#             points = np.column_stack([lt, ll])
# 
#             # get interferogram full grid
#             lats = grid_ll.lat.values
#             lons = grid_ll.lon.values
# 
#             # calculate trans grid subset extent
#             ltmin, ltmax = np.nanmin(lt), np.nanmax(lt)
#             llmin, llmax = np.nanmin(ll), np.nanmax(ll)
#             # and spacing
#             dlat = np.diff(lats)[0]
#             dlon = np.diff(lons)[0]
# 
#             # select required interferogram grid subset
#             lats = lats[(lats>ltmin-dlat)&(lats<ltmax+dlat)]
#             lons = lons[(lons>llmin-dlon)&(lons<llmax+dlon)]
# 
#             # for cropped interferogram we can have no valid pixels for the processing
#             if lats.size == 0 or lons.size == 0:
#                 return np.nan * np.zeros((trans_inv.y.size, trans.x.size), dtype=np.float32)
# 
#             values = grid_ll.sel(lat=lats, lon=lons).values.astype(np.float64)
# 
#             # perform interpolation
#             interp = RegularGridInterpolator((lats, lons), values, method='nearest', bounds_error=False)
#             grid_ra = interp(points).reshape(trans_inv.y.size, trans_inv.x.size).astype(np.float32)
# 
#             return grid_ra
# 
#         # get transform table
#         trans_inv = self.get_intf_ll2ra()[['lt', 'll']]
#         ys = trans_inv.y
#         xs = trans_inv.x
#         # define processing blocks
#         chunks = ys.size / chunksize
#         xs_blocks = np.array_split(xs, np.arange(0, xs.size, chunksize)[1:])
# 
#         grids_ra = []
#         # unify input grid(s) to stack
#         for grid_ll in grids if len(grids.dims) == 3 else [grids]:
#             # per-block processing
#             blocks  = []
#             for xs_block in xs_blocks:
#                 block = dask.array.from_delayed(intf_block_inv(trans_inv.sel(x=xs_block), grid_ll),
#                                                 shape=(ys.size, xs_block.size), dtype=np.float32)
#                 blocks.append(block)
#             #grid_ll = dask.array.block(blocks)
#             # set the output grid and drop the fake dimension if needed
#             grid_ra = xr.DataArray(dask.array.block(blocks), coords=trans_inv.coords).rename(grids.name)
#             grids_ra.append(grid_ra)
# 
#         if len(grids.dims) == 2:
#             # drop the fake dimension
#             coords = trans_inv.coords
#             out = xr.DataArray(grids_ra[0], coords=coords).rename(grids.name)
#         else:
#             # find stack dim
#             stack_dim = grids.dims[0]
#             coords = {stack_dim: grids[stack_dim], 'y': trans_inv.y, 'x': trans_inv.x}
#             out = xr.DataArray(dask.array.stack(grids_ra), coords=coords).rename(grids.name)
# 
#         # append source grid coordinates excluding removed lat, lon ones
#         for (k,v) in grids.coords.items():
#             if k not in ['lat','lon']:
#                 out[k] = v
#         return out.chunk(chunksize)
