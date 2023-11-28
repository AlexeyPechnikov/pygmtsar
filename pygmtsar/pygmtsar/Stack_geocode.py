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

    def compute_geocode(self, coarsen=60.):
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
        # suppress Dask warning "RuntimeWarning: All-NaN slice encountered"
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')

        self.compute_trans(coarsen=coarsen)
        self.compute_trans_inv(coarsen=coarsen)
        self.compute_satellite_look_vector()

    # coarsen=1:
    # nearest: coords [array([596.97436523]), array([16976.93164062])]
    # linear:  coords [array([597.10408125]), array([16977.41447869])]
    # coarsen=2:
    # nearest: coords [array([596.4810791]), array([16977.52734375])]
    # linear:  coords [array([597.10471665]), array([16977.3737493])]
    # coarsen=4:
    # nearest: coords [array([596.42352295]), array([16978.65625])]
    # linear:  coords [array([597.1080563]), array([16977.35608873])]
    def geocode(self, geometry, z_offset=None):
        """
        Inverse geocode input geodataframe with 2D or 3D points. 
        
        Examples:
        ---------
        sbas.geocode(AOI.assign(col=0))
        sbas.geocode(gpd.GeoDataFrame(geometry=AOI.centroid))
        sbas.geocode(AOI.centroid)
        sbas.geocode(AOI.geometry.item())
        sbas.geocode([AOI.geometry.item()])
        sbas.geocode(geometry)
        """
        import numpy as np
        import geopandas as gpd
        import xarray as xr
        from shapely.geometry import LineString, Point

        if isinstance(geometry, (gpd.GeoDataFrame, gpd.GeoSeries)):
            geometries = geometry.geometry
        elif isinstance(geometry, (list, tuple, np.ndarray)):
            geometries = geometry
        else:
            geometries = [geometry]

        dem = self.get_dem()
        prm = self.PRM()
    
        geoms = []
        for geom in geometries:
            #print (geom)
            coords = np.asarray(geom.coords[:])
            #print (len(coords))
            ele = dem.interp(lat=xr.DataArray(coords[:,1]),
                       lon=xr.DataArray(coords[:,0]), method='linear').compute()
            if z_offset is None:
                z = coords[:,2] if coords.shape[1]==3 else 0
            else:
                z = z_offset
            points_ll = np.column_stack([coords[:,0], coords[:,1], ele.values + z])
            #print (points_ll)
            points_ra = prm.SAT_llt2rat(points_ll)
            points_ra = points_ra[:,:2] if len(points_ll)>1 else [points_ra[:2]]
            #print (points_ra)
            geom_type = geom.type
            if geom_type == 'LineString':
                geom = LineString(points_ra)
            elif geom_type == 'Point':
                geom = Point(points_ra)
            #print (geom)
            geoms.append(geom)

        # preserve original dataframe columns if exists
        if isinstance(geometry, gpd.GeoDataFrame):
            return gpd.GeoDataFrame(geometry, geometry=geoms, crs=3857)
        elif isinstance(geometry, gpd.GeoSeries):
            return gpd.GeoSeries(geoms, crs=3857)
        elif isinstance(geometry, (list, tuple, np.ndarray)):
            return geoms
        else:
            return geoms[0]

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
            trans = trans.sel(lat=trans.lat[step_y//2::step_y], lon=trans.lon[step_x//2::step_x])
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
    def ll2ra(self, data, autoscale=True):
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

        # analyse grid and transform matrix spacing
        grid_dlat = np.diff(data.lat)[0]
        grid_dlon = np.diff(data.lon)[0]
        trans_inv_dlat = np.diff(trans_inv.lat)[0]
        trans_inv_dlon = np.diff(trans_inv.lon)[0]

        # define transform spacing in radar coordinates
        step_lat = int(np.round(grid_dlat / trans_inv_dlat))
        step_lon = int(np.round(grid_dlon / trans_inv_dlon))
        #print ('step_lat', step_lat, 'step_lon', step_lon)
        assert step_lat>=1 and step_lon>=1, f'Transforming grid spacing (grid_dlat, grid_dlon) is smaller \
                                          than transform matrix spacing (trans_inv_dlat, trans_inv_dlon), \
                                          call Stack.geocode() with less coarsing'
        # decimate the full trans_inv grid to the required spacing
        if autoscale and (step_lat>1 or step_lon>1):
            # define the equally spacing geographic coordinates grid
            trans_inv = trans_inv.sel(y=trans_inv.y[step_lat//2::step_lat], x=trans_inv.x[step_lon//2::step_lon])
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
