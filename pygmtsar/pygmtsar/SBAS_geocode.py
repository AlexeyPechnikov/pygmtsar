# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_sbas import SBAS_sbas
from .tqdm_dask import tqdm_dask

class SBAS_geocode(SBAS_sbas):


    # coarsen=1:
    # nearest: coords [array([596.97436523]), array([16976.93164062])]
    # linear:  coords [array([597.10408125]), array([16977.41447869])]
    # coarsen=2:
    # nearest: coords [array([596.4810791]), array([16977.52734375])]
    # linear:  coords [array([597.10471665]), array([16977.3737493])]
    # coarsen=4:
    # nearest: coords [array([596.42352295]), array([16978.65625])]
    # linear:  coords [array([597.1080563]), array([16977.35608873])]
    def ll2ra(self, data):
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
            ele = interp([geom.y, geom.x])[0] + (geom.z if geom.has_z else 0)
            points_ll.append([geom.x, geom.y, ele])
            del interp
        points_ra = prm.SAT_llt2rat(points_ll)[:,:2]
        #print ('points_ra', points_ra)
        # set fake CRS to differ from WGS84 coordinates
        return gpd.GeoDataFrame(data, geometry=[Point(point_ra) for point_ra in points_ra], crs=3857)

    def geocode_parallel(self, pairs=None, coarsen=None, chunksize=None):
        """
        Perform parallel geocoding of the interferograms.

        Parameters
        ----------
        pairs : list or None, optional
            List of interferogram pairs to process. If None, all available pairs will be processed.

        Notes
        -----
        This method performs parallel geocoding of the interferograms. It builds the necessary geocoding matrices
        to apply the geocoding transformation to the interferogram grids. The geocoding involves converting the
        radar coordinates to geographic coordinates and vice versa.
        """

        # build trans_dat, trans_ dat_inv and topo_ra grids for merged subswaths
        # for single-digit subswath the grids already created for interferogram processing
        if len(str(self.get_subswath())) > 1 or coarsen is not None:
            self.topo_ra_parallel(coarsen=coarsen)

        # build geographic coordinates transformation matrix for landmask and other grids
        self.intf_ll2ra_matrix_parallel(pairs=pairs, chunksize=chunksize)
        # build radar coordinates transformation matrix for the interferograms grid stack        
        self.intf_ra2ll_matrix_parallel(pairs=pairs, chunksize=chunksize)

##########################################################################################
# ra2ll
##########################################################################################
    def intf_ra2ll_matrix_parallel(self, pairs=None, chunksize=None, interactive=False):
        """
        Perform parallel computation of the radar-to-geographic coordinate transformation matrix.

        Parameters
        ----------
        pairs : list or None, optional
            List of interferogram pairs to process. If None, all available pairs will be processed.
        interactive : bool, optional
            Flag indicating whether to return the matrix without saving to a file.
        """
        import xarray as xr
        import numpy as np
        import dask
        import os

        if chunksize is None:
            chunksize = self.chunksize

        pairs = self.pairs(pairs)

        # find any one interferogram to define the grid
        intf = self.open_grids(pairs[:1], 'phasefilt')[0].astype(bool)
        dy = np.diff(intf.y)[0]
        dx = np.diff(intf.x)[0]

        # get transform table
        trans = self.get_trans_dat()
        trans_dy = np.diff(trans.y)[0]
        trans_dx = np.diff(trans.x)[0]

        # define transform spacing in radar coordinates
        step_y = int(np.round(dy / trans_dy))
        step_x = int(np.round(dx / trans_dx))
        #print ('step_y', step_y, 'step_x', step_x)

        assert step_y>=1 and step_x>=1, 'Interferogram grid spacing is smaller than transform grid spacing; \
                                          call SBAS.geocode_parallel() with less coarsing'

        # define the equally spacing geographic coordinates grid
        lats = trans.lat[::step_y]
        lons = trans.lon[::step_x]

        # decimate the full trans grid to the required spacing
        trans = trans.sel(lat=lats, lon=lons)[['azi', 'rng']]
        # define interferogram radar coordinates grid
        trans['y'] = xr.DataArray(intf.y.values, dims='y')
        trans['x'] = xr.DataArray(intf.x.values, dims='x')

        if interactive:
            return trans

        # save to NetCDF file
        filename = self.get_filenames(None, None, 'intf_ra2ll')
        #print ('filename', filename)
        # to resolve NetCDF rewriting error
        if os.path.exists(filename):
            os.remove(filename)
        encoding = {var: self.compression(trans[var].shape, chunksize=chunksize) for var in trans.data_vars}
        handler = trans.to_netcdf(filename,
                                  encoding=encoding,
                                  engine=self.engine,
                                  compute=False)
        tqdm_dask(dask.persist(handler), desc='Build ra2ll Transform')
        # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
        import gc; gc.collect()

    def get_intf_ra2ll(self, subswath=None, chunksize=None):
        """
        Get the radar-to-geographic coordinate transformation matrix.

        Parameters
        ----------
        chunksize : int or dict, optional
            Chunk size for dask arrays. If not provided, the default chunk size is used.

        Returns
        -------
        xarray.DataArray
            The radar-to-geographic coordinate transformation matrix.

        Notes
        -----
        This method retrieves the radar-to-geographic coordinate transformation matrix (intf_ra2ll) stored in the
        NetCDF grid. The matrix is useful for inverse geocoding, converting interferogram grids from radar
        coordinates to geographic coordinates.
        """
        import xarray as xr

        subswath = self.get_subswath(subswath)
        filename = self.get_filenames(subswath, None, 'intf_ra2ll')
        trans = xr.open_dataset(filename, engine=self.engine, chunks=self.chunksize)
        return trans

    def intf_ra2ll(self, grid, chunksize=None):
        """
        Faster geocoding for interferogram grid. It uses usually decimated transformation table vs full table in ra2ll.
        """
        trans = self.get_intf_ra2ll()
        return self.ra2ll(grid, trans=trans, chunksize=chunksize)

    # intf_ra2ll
    def ra2ll(self, grid, trans=None, chunksize=None):
        """
        Perform geocoding from radar to geographic coordinates.

        Parameters
        ----------
        grid : xarray.DataArray
            Grid(s) representing the interferogram(s) in radar coordinates.
        trans : xarray.DataArray
            Geocoding transform matrix in radar coordinates.
        chunksize : int or dict, optional
            Chunk size for dask arrays. If not provided, the default chunk size is used.

        Returns
        -------
        xarray.DataArray
            The inverse geocoded grid(s) in geographic coordinates.

        Examples
        --------
        Geocode 3D unwrapped phase grid stack:
        unwraps_ll = sbas.intf_ra2ll(sbas.open_grids(pairs, 'unwrap'))
        # or use "geocode" option for open_grids() instead:
        unwraps_ll = sbas.open_grids(pairs, 'unwrap', geocode=True)
        """
        import dask
        import xarray as xr
        import numpy as np

        # helper check
        if not 'y' in grid.dims or not 'x' in grid.dims:
            print ('NOTE: the input grid is not in radar coordinates, miss geocoding')
            return grid

        if chunksize is None:
            chunksize = self.chunksize

        if trans is None:
            # get complete transform table
            trans = self.get_trans_dat()

        @dask.delayed
        def intf_block(lats_block, lons_block, grid_ra):
            from scipy.interpolate import RegularGridInterpolator

            trans_block = trans.sel(lat=lats_block, lon=lons_block)
            # check if the data block exists
            if not (trans_block.lat.size>0 and trans_block.lon.size>0):
                return np.nan * np.zeros((lats_block.size, lons_block.size), dtype=np.float32)

            # use trans table subset
            y = trans_block.azi.values.ravel()
            x = trans_block.rng.values.ravel()
            points = np.column_stack([y, x])

            # get interferogram full grid
            ys = grid_ra.y.values
            xs = grid_ra.x.values

            # calculate trans grid subset extent
            ymin, ymax = np.nanmin(y), np.nanmax(y)
            xmin, xmax = np.nanmin(x), np.nanmax(x)
            # and spacing
            dy = np.diff(ys)[0]
            dx = np.diff(xs)[0]

            # select required interferogram grid subset
            ys = ys[(ys>ymin-dy)&(ys<ymax+dy)]
            xs = xs[(xs>xmin-dx)&(xs<xmax+dx)]

            # for cropped interferogram we can have no valid pixels for the processing
            if ys.size == 0 or xs.size == 0:
                return np.nan * np.zeros((lats_block.size, lons_block.size), dtype=np.float32)

            values = grid_ra.sel(y=ys, x=xs).values.astype(np.float64)

            # perform interpolation
            interp = RegularGridInterpolator((ys, xs), values, method='nearest', bounds_error=False)
            grid_ll = interp(points).reshape(lats_block.size, lons_block.size).astype(np.float32)

            return grid_ll

        # analyse grid and transform matrix spacing
        grid_dy = np.diff(grid.y)[0]
        grid_dx = np.diff(grid.x)[0]
        trans_dy = np.diff(trans.y)[0]
        trans_dx = np.diff(trans.x)[0]

        # define transform spacing in radar coordinates
        step_y = int(np.round(grid_dy / trans_dy))
        step_x = int(np.round(grid_dx / trans_dx))
        #print ('step_y', step_y, 'step_x', step_x)

        assert step_y>=1 and step_x>=1, f'Transforming grid spacing (grid_dy, grid_dx) is smaller \
                                          than transform matrix spacing (trans_dy, trans_dx), \
                                          call SBAS.topo_ra_parallel() or SBAS.geocode_parallel() with less coarsing'
        # decimate the full trans grid to the required spacing
        if step_y>1 or step_x>1:
            # define the equally spacing geographic coordinates grid
            trans = trans.sel(lat=trans.lat[::step_y], lon=trans.lon[::step_x])
        # select required variables only
        trans = trans[['azi', 'rng']]

        # define the equally spacing geographic coordinates grid
        lats = trans.lat[::step_y]
        lons = trans.lon[::step_x]

        # specify output coordinates
        lats = trans.lat
        lons = trans.lon
        # split to equal chunks and rest
        lats_blocks = np.array_split(lats, np.arange(0, lats.size, chunksize)[1:])
        lons_blocks = np.array_split(lons, np.arange(0, lons.size, chunksize)[1:])

        # find stack dim
        stackvar = grid.dims[0] if len(grid.dims) == 3 else 'stack'
        #print ('stackvar', stackvar)

        stack = []
        for stackval in grid[stackvar].values if len(grid.dims) == 3 else [None]:
            block_grid = grid.sel({stackvar: stackval}) if stackval is not None else grid
            # per-block processing
            blocks_total  = []
            for lats_block in lats_blocks:
                blocks = []
                for lons_block in lons_blocks:
                    block = dask.array.from_delayed(intf_block(lats_block, lons_block, block_grid),
                                                    shape=(lats_block.size, lons_block.size), dtype=np.float32)
                    blocks.append(block)
                    del block
                blocks_total.append(blocks)
                del blocks
            dask_block = dask.array.block(blocks_total)
            if len(grid.dims) == 3:
                coords = {stackvar: [stackval], 'lat': trans.coords['lat'], 'lon': trans.coords['lon']}
                da = xr.DataArray(dask_block[None, :], coords=coords)
            else:
                coords = {'lat': trans.coords['lat'], 'lon': trans.coords['lon']}
                da = xr.DataArray(dask_block, coords=coords)
            stack.append(da)
            del blocks_total, dask_block

        # wrap lazy Dask array to Xarray dataarray
        if len(grid.dims) == 2:
            out = stack[0]
        else:
            out = xr.concat(stack, dim=stackvar)
        del stack

        # append source grid coordinates excluding removed y, x ones
        for (k,v) in grid.coords.items():
            if k not in ['y','x']:
                out[k] = v
        return out.rename(grid.name)
    
##########################################################################################
# ll2ra
##########################################################################################
    def intf_ll2ra_matrix_parallel(self, pairs=None, chunksize=None, interactive=False):
        """
        Perform parallel computation of the geographic-to-radar coordinate transformation matrix.

        Parameters
        ----------
        pairs : list or None, optional
            List of interferogram pairs to process. If None, all available pairs will be processed.
        interactive : bool, optional
            Flag indicating whether to return the matrix without saving to a file.
        """
        import xarray as xr
        import numpy as np
        import dask
        import os

        if chunksize is None:
            chunksize = self.chunksize

        pairs = self.pairs(pairs)

        @dask.delayed
        def intf_block_inv(intf_block):
            from scipy.interpolate import RegularGridInterpolator

            # define input grid subset
            trans_inv_dy = np.diff(trans_inv.y)[0]
            trans_inv_dx = np.diff(trans_inv.x)[0]
            ymin, ymax = np.min(intf_block.y)-trans_inv_dy, np.max(intf_block.y)+trans_inv_dy
            xmin, xmax = np.min(intf_block.x)-trans_inv_dx, np.max(intf_block.x)+trans_inv_dx
            trans_block = trans_inv.sel(y=slice(ymin, ymax), x=slice(xmin,xmax))
            #print ('trans_block', trans_block)

            # for cropped interferogram we can have no valid pixels for the processing
            if trans_block.y.size == 0 or trans_block.x.size == 0:
                fake_block = np.nan * np.zeros((intf.y.size, intf.x.size), dtype=np.float32)
                return np.asarray([fake_block, fake_block])

            # define input grid
            trans_grid = (trans_block.y, trans_block.x)

            # define output grid
            ys, xs = np.meshgrid(intf_block.y, intf_block.x)
            points = np.column_stack([ys.ravel(), xs.ravel()])

            grids = []
            for var in ['lt', 'll']:
                # np.float64 required for the interpolator
                block = trans_block[var].compute(n_workers=1).data.astype(np.float64)
                # perform interpolation
                interp = RegularGridInterpolator(trans_grid, block, method='linear', bounds_error=False)
                grid = interp(points).reshape((intf_block.y.size, intf_block.x.size))
                grids.append(grid)
                del block, interp, grid

            del trans_block, trans_grid, ys, xs, points
            return np.asarray(grids)

        # find any one interferogram to define the grid
        intf = self.open_grids(pairs[:1], 'phasefilt')[0].astype(bool)
        #print ('intf', intf)
        blocks = chunksize**2 // intf.y.size
        xs_blocks = np.array_split(intf.x, np.arange(0, intf.x.size, blocks if blocks>=1 else 1)[1:])
        #print ('xs_blocks', xs_blocks)

        # get transform table
        trans_inv = self.get_trans_dat_inv()[['lt', 'll']]
        #print ('trans_inv', trans_inv)

        # per-block processing
        block_lts  = []
        block_lls  = []
        for xs_block in xs_blocks:
            block_lt, block_ll = dask.array.from_delayed(intf_block_inv(intf.sel(x=xs_block)),
                                            shape=(2, intf.y.size, xs_block.size), dtype=np.float32)
            block_lts.append(block_lt)
            block_lls.append(block_ll)
            del block_lt, block_ll
        # set the output grid and drop the fake dimension if needed
        lt = xr.DataArray(dask.array.block(block_lts), coords=intf.coords)
        ll = xr.DataArray(dask.array.block(block_lls), coords=intf.coords)
        trans = xr.Dataset({'lt': lt, 'll': ll}).rename({'y': 'a', 'x': 'r'})
    
        if interactive:
            return trans

        # save to NetCDF file
        filename = self.get_filenames(None, None, 'intf_ll2ra')
        #print ('filename', filename)
        # to resolve NetCDF rewriting error
        if os.path.exists(filename):
            os.remove(filename)
        encoding = {var: self.compression(trans[var].shape, chunksize=chunksize) for var in trans.data_vars}
        handler = trans.to_netcdf(filename,
                                  encoding=encoding,
                                  engine=self.engine,
                                  compute=False)
        tqdm_dask(dask.persist(handler), desc='Build ll2ra Transform')
        # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
        import gc; gc.collect()

    def get_intf_ll2ra(self, subswath=None, chunksize=None):
        """
        Get the geographic-to-radar coordinate transformation matrix.

        Parameters
        ----------
        chunksize : int or dict, optional
            Chunk size for dask arrays. If not provided, the default chunk size is used.

        Returns
        -------
        xarray.DataArray
            The radar-to-geographic coordinate transformation matrix.

        Notes
        -----
        This method retrieves the geographic-to-radar coordinate transformation matrix (intf_ll2ra) stored in the
        NetCDF grid. The matrix is useful for direct geocoding, converting geographic coordinate grids to
        radar coordinates interferogram grid.
        """
        import xarray as xr

        subswath = self.get_subswath(subswath)
        filename = self.get_filenames(subswath, None, 'intf_ll2ra')
        trans_inv = xr.open_dataset(filename, engine=self.engine, chunks=self.chunksize)
        if 'a' in trans_inv.dims and 'r' in trans_inv.dims:
            return trans_inv.rename({'a': 'y', 'r': 'x'})
        else:
            return trans_inv

    def intf_ll2ra(self, grids, chunksize=None):
        """
        Perform inverse geocoding from geographic to radar coordinates.

        Parameters
        ----------
        grids : xarray.DataArray
            Grid(s) representing the interferogram(s) in radar coordinates.
        chunksize : int or dict, optional
            Chunk size for dask arrays. If not provided, the default chunk size is used.

        Returns
        -------
        xarray.DataArray
            The inverse geocoded grid(s) in geographic coordinates.

        Examples
        --------
        Inverse geocode 3D unwrapped phase grids stack:
        unwraps_ll = sbas.open_grids(pairs, 'unwrap', geocode=True)
        unwraps = sbas.intf_ll2ra(unwraps_ll)
        """
        import dask
        import xarray as xr
        import numpy as np

        # helper check
        if not 'lat' in grids.dims or not 'lon' in grids.dims:
            print ('NOTE: the input grid is not in geographic coordinates, miss inverse geocoding')
            return grids

        if chunksize is None:
            chunksize = self.chunksize

        @dask.delayed
        def intf_block_inv(trans_inv, grid_ll):
            from scipy.interpolate import RegularGridInterpolator

            # use trans table subset
            lt = trans_inv.lt.values.ravel()
            ll = trans_inv.ll.values.ravel()
            points = np.column_stack([lt, ll])

            # get interferogram full grid
            lats = grid_ll.lat.values
            lons = grid_ll.lon.values

            # calculate trans grid subset extent
            ltmin, ltmax = np.nanmin(lt), np.nanmax(lt)
            llmin, llmax = np.nanmin(ll), np.nanmax(ll)
            # and spacing
            dlat = np.diff(lats)[0]
            dlon = np.diff(lons)[0]

            # select required interferogram grid subset
            lats = lats[(lats>ltmin-dlat)&(lats<ltmax+dlat)]
            lons = lons[(lons>llmin-dlon)&(lons<llmax+dlon)]

            # for cropped interferogram we can have no valid pixels for the processing
            if lats.size == 0 or lons.size == 0:
                return np.nan * np.zeros((trans_inv.y.size, trans.x.size), dtype=np.float32)

            values = grid_ll.sel(lat=lats, lon=lons).values.astype(np.float64)

            # perform interpolation
            interp = RegularGridInterpolator((lats, lons), values, method='nearest', bounds_error=False)
            grid_ra = interp(points).reshape(trans_inv.y.size, trans_inv.x.size).astype(np.float32)

            return grid_ra

        # get transform table
        trans_inv = self.get_intf_ll2ra()[['lt', 'll']]
        ys = trans_inv.y
        xs = trans_inv.x
        # define processing blocks
        chunks = ys.size / chunksize
        xs_blocks = np.array_split(xs, np.arange(0, xs.size, chunksize)[1:])

        grids_ra = []
        # unify input grid(s) to stack
        for grid_ll in grids if len(grids.dims) == 3 else [grids]:
            # per-block processing
            blocks  = []
            for xs_block in xs_blocks:
                block = dask.array.from_delayed(intf_block_inv(trans_inv.sel(x=xs_block), grid_ll),
                                                shape=(ys.size, xs_block.size), dtype=np.float32)
                blocks.append(block)
            #grid_ll = dask.array.block(blocks)
            # set the output grid and drop the fake dimension if needed
            grid_ra = xr.DataArray(dask.array.block(blocks), coords=trans_inv.coords).rename(grids.name)
            grids_ra.append(grid_ra)

        if len(grids.dims) == 2:
            # drop the fake dimension
            coords = trans_inv.coords
            out = xr.DataArray(grids_ra[0], coords=coords).rename(grids.name)
        else:
            # find stack dim
            stack_dim = grids.dims[0]
            coords = {stack_dim: grids[stack_dim], 'y': trans_inv.y, 'x': trans_inv.x}
            out = xr.DataArray(dask.array.stack(grids_ra), coords=coords).rename(grids.name)

        # append source grid coordinates excluding removed lat, lon ones
        for (k,v) in grids.coords.items():
            if k not in ['lat','lon']:
                out[k] = v
        return out.chunk(chunksize)
