# ----------------------------------------------------------------------------
# PyGMTSAR
#.
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
#
# Copyright (c) 2024, Alexey Pechnikov
#
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_ps import Stack_ps
from .tqdm_dask import tqdm_dask

class Stack_export(Stack_ps):

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
        stack.as_geo(grid).rio.clip([geometry])

        Notes
        -----
        This method adds geospatial attributes (CRS and spatial dimensions) to the input grid,
        allowing raster operations using the RioXarray library. If the input grid is already
        in geographic coordinates, the CRS is set to EPSG 4326 with spatial dimensions 'lat' and 'lon'.
        Otherwise, if the input grid is in radar coordinates, a fake metric coordinate system is used
        with EPSG 3857 and spatial dimensions 'y' and 'x'. The method relies on the availability of the
        'rioxarray' module.
        """
        import rioxarray
        import sys
        #assert 'rioxarray' in sys.modules, 'rioxarray module is not found'
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
    def as_vtk(dataset):
        """
        Map a 2D image onto a 3D topography.

        Parameters:
            dataset: A file name or Xarray Dataset to use for mapping the image.

        Returns:
            A vtkStructuredGrid representing the image mapped onto the topography.
        """
        from vtk import vtkPoints, vtkStructuredGrid, vtkThreshold, vtkDataObject, \
            VTK_FLOAT, VTK_UNSIGNED_CHAR, vtkStringArray, vtkFloatArray, vtkIntArray
        from vtk.util import numpy_support as vn
        import numpy as np
        import xarray as xr

        assert isinstance(dataset, xr.Dataset), 'ERROR: Expected "dataset" argument as Xarray Dataset'

        # fill NODATA by NAN
        for data_var in dataset.data_vars:
            #if dataset[data_var].values.dtype in [np.dtype('float16'),np.dtype('float32'),np.dtype('float64')]:
            if np.issubdtype(dataset[data_var].dtype, np.floating):
                dataset[data_var].values = dataset[data_var].values.astype('float32')
                if '_FillValue' in dataset[data_var] and not np.isnan(dataset[data_var]._FillValue):
                    dataset[data_var].values[dataset[data_var].values == dataset[data_var]._FillValue] = np.nan

        xs = dataset.x.values
        ys = dataset.y.values
        (yy,xx) = np.meshgrid(ys, xs)
        if 'z' in dataset:
            zs = dataset.z.values
        else:
            # set z coordinate to zero
            zs = np.zeros_like(xx)

        # create raster mask by geometry and for NaNs
        vtk_points = vtkPoints()
        points = np.column_stack((xx.ravel('F'), yy.ravel('F'), zs.ravel('C')))
        vtk_points.SetData(vn.numpy_to_vtk(points, deep=True))

        sgrid = vtkStructuredGrid()
        sgrid.SetDimensions(len(xs), len(ys), 1)
        sgrid.SetPoints(vtk_points)

        for data_var in dataset.data_vars:
            if 'z' == data_var:
                # variable is already included into coordinates
                continue
            da = dataset[data_var]
            dims = da.dims
            #print (data_var, dims, da.dtype)
            if 'band'in dims:
                bands = da.band.shape[0]
                #print ('bands', bands)
                if bands in [3,4]:
                    # RGB or RGBA, select 3 bands only
                    array = vn.numpy_to_vtk(da[:3].round().astype(np.uint8).values.reshape(3,-1).T, deep=True, array_type=VTK_UNSIGNED_CHAR)
                    array.SetName(da.name)
                    sgrid.GetPointData().AddArray(array)
                elif bands == 1:
                    array = vn.numpy_to_vtk(da.values.reshape(1,-1).T, deep=True, array_type=VTK_FLOAT)
                    array.SetName(da.name)
                    sgrid.GetPointData().AddArray(array)
                else:
                    print (f'ERROR: Unsupported bands count (should be 1,3 or 4) {bands} for variable {data_var}')
                    return
            else:
                array = vn.numpy_to_vtk(da.values.ravel(), deep=True, array_type=VTK_FLOAT)
                array.SetName(da.name)
                sgrid.GetPointData().AddArray(array)

        for coord in dataset.coords:
            #print (coord, dataset[coord].dims, len(dataset[coord].dims))
            if len(dataset[coord].dims) > 0:
                # real coordinate, ignore it
                continue
            #print (coord, dataset[coord].dtype)
            #print ('np.datetime64', np.issubdtype(dataset[coord].dtype, np.datetime64))
            #print ('np.int64', np.issubdtype(dataset[coord].dtype, np.int64))
            #print ('np.float64', np.issubdtype(dataset[coord].dtype, np.float64))
            if np.issubdtype(dataset[coord].dtype, np.datetime64):
                data_array = vtkStringArray()
                data_value = str(dataset[coord].dt.date.values)
            elif np.issubdtype(dataset[coord].dtype, str):
                data_array = vtkStringArray()
                data_value = str(dataset[coord].values)
            elif np.issubdtype(dataset[coord].dtype, np.int64):
                data_array = vtkIntArray()
                data_value = dataset[coord].values.astype(np.int64)
            elif np.issubdtype(dataset[coord].dtype, np.float64):
                data_array = vtkFloatArray
                data_value = dataset[coord].values.astype(np.float64)
            else:
                print(f'NOTE: unsupported attribute {coord} datatype {dataset[coord].dtype}, miss it')
                continue
            data_array.SetName(coord)
            data_array.InsertNextValue(data_value)
            sgrid.GetFieldData().AddArray(data_array)

        # required for 3D interactive rendering on Google Colab
        return sgrid

    def export_geotiff(self, data, name, caption='Exporting WGS84 GeoTIFF(s)'):
        """
        Export single GeoTIFF file "velocity.tif":
        sbas.export_geotiff(velocity, 'velocity')
        
        Export pair-based stack of GeoTIFF files like corr.2024-01-01_2024-01-02.tif, ...:
        sbas.export_geotiff(corr, 'corr')
        
        Export date-based stack of GeoTIFF files (ike disp.2024-01-01.tif, ...:
        sbas.export_geotiff(disp, 'disp')

        The exported GeoTIFF files can be converted to KMZ for Google Earth Engine using GDAL tools:
        gdalwarp -of KMLSUPEROVERLAY -co FORMAT=PNG velocity.tif velocity.kmz
        """
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm

        assert isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object'

        # determine if data has a stack dimension and what it is
        stackvar = data.dims[0] if len(data.dims) == 3 else None

        # prepare the progress bar
        with tqdm(desc=caption, total=len(data[stackvar]) if stackvar is not None else 1) as pbar:
            for grid in data.transpose(stackvar, ...) if stackvar is not None else [data]:
                # convert the stack variable value to a string suitable for filenames
                if stackvar is not None and np.issubdtype(grid[stackvar].dtype, np.datetime64):
                    stackval = grid[stackvar].dt.date.item()
                elif stackvar is not None:
                    stackval = grid[stackvar].astype(str).item().replace(' ', '_')
                else:
                    stackval = ''
                #print ('stackval', stackval)
                filename = f'{name}.{stackval}.tif' if stackvar is not None else f'{name}.tif'
                #print ('filename', filename)
                # convert the data to geographic coordinates if necessary and export to GeoTIFF
                self.as_geo(self.ra2ll(grid) if not self.is_geo(grid) else grid).rio.to_raster(filename)
                pbar.update(1)

    def export_geojson(self, data, name, caption='Exporting WGS84 GeoJSON', pivotal=True, digits=1):
        """
        Export GeoJSON file "velocity.geojson":
        sbas.export_geojson(velocity, 'velocity')
        """
        import xarray as xr
        import geopandas as gpd
        import numpy as np
        from tqdm.auto import tqdm
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        from dask.distributed import utils_perf
        utils_perf.disable_gc_diagnosis()

        assert isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object'

        # determine if data has a stack dimension and what it is
        stackvar = data.dims[0] if len(data.dims) == 3 else None

        # convert the data to geographic coordinates if necessary
        if not self.is_geo(data):
            grid = self.ra2ll(data)
        else:
            grid = data

        # split to equal chunks and rest
        lats_blocks = np.array_split(np.arange(grid.lat.size), np.arange(0, grid.lat.size, self.netcdf_chunksize//2)[1:])
        lons_blocks = np.array_split(np.arange(grid.lon.size), np.arange(0, grid.lon.size, self.netcdf_chunksize//2)[1:])

        empty = True
        # prepare the progress bar
        with tqdm(desc=caption, total=len(lats_blocks)*len(lons_blocks)) as pbar:
            with open(f'{name}.geojson', 'w') as f:
                # GeoJSON header
                f.write('{"type": "FeatureCollection", "features": [\n')
                for lats_block in lats_blocks:
                    for lons_block in lons_blocks:
                        block = grid.isel(lat=lats_block, lon=lons_block).compute()
                        df = block.to_dataframe().dropna().reset_index()
                        del block
                        df[data.name] = df[data.name].apply(lambda x: float(f"{x:.{digits}f}"))
                        for col in df.columns:
                            if np.issubdtype(df[col].dtype, np.datetime64):
                                df[col] = df[col].dt.date.astype(str)
                        if stackvar is not None and np.issubdtype(df[stackvar].dtype, np.datetime64):
                            df[stackvar] = df[stackvar].dt.date.astype(str)
                        if stackvar is not None and pivotal:
                            df = df.pivot_table(index=['lat', 'lon'], columns=stackvar,
                                                values=data.name, fill_value=np.nan).reset_index()
                        # convert to geodataframe with value and geometry columns
                        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon.round(6), df.lat.round(6)))
                        del df, gdf['lat'], gdf['lon']
                        if len(gdf):
                            # crop GeoJSON header and footer and split lines
                            chunk_str = gdf.to_json(drop_id=True)[43:-2].replace('}}, ', '}},\n')
                            f.write(chunk_str if empty else ',' + chunk_str)
                            empty = False
                            del chunk_str
                        del gdf
                        pbar.update(1)
                # GeoJSON footer
                f.write(']}')
        del grid

    def export_csv(self, data, name, caption='Exporting WGS84 CSV', digits=1, delimiter=','):
        """
        Export CSV file "velocity.csv":
        sbas.export_csv(velocity, 'velocity')
        """
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        from dask.distributed import utils_perf
        utils_perf.disable_gc_diagnosis()
    
        assert isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object'
    
        # determine if data has a stack dimension and what it is
        stackvar = data.dims[0] if len(data.dims) == 3 else None
    
        # convert the data to geographic coordinates if necessary
        if not self.is_geo(data):
            grid = self.ra2ll(data)
        else:
            grid = data
    
        # split to equal chunks and rest
        lats_blocks = np.array_split(np.arange(grid.lat.size), np.arange(0, grid.lat.size, self.netcdf_chunksize//2)[1:])
        lons_blocks = np.array_split(np.arange(grid.lon.size), np.arange(0, grid.lon.size, self.netcdf_chunksize//2)[1:])
    
        # prepare the progress bar
        with tqdm(desc=caption, total=len(lats_blocks)*len(lons_blocks)) as pbar:
            with open(f'{name}.csv', 'w') as f:
                # CSV header
                f.write(delimiter.join(filter(None, [stackvar, 'lon', 'lat', data.name])) + '\n')
                for lats_block in lats_blocks:
                    for lons_block in lons_blocks:
                        block = grid.isel(lat=lats_block, lon=lons_block).compute()
                        block_val = block.round(digits).values
                        block_lat = block.lat.round(6).values
                        block_lon = block.lon.round(6).values
                        if stackvar is not None:
                            stackvals = block[stackvar]
                            if np.issubdtype(stackvals.dtype, np.datetime64):
                                stackvals = stackvals.dt.date.astype(str)                           
                            stackvals, lats, lons = np.meshgrid(stackvals, block_lat, block_lon, indexing='ij')
                            block_csv = np.column_stack((stackvals.ravel(), lons.ravel(), lats.ravel(), block_val.ravel()))
                            del stackvals, lats, lons
                        else:
                            lats, lons = np.meshgrid(block_lat, block_lon, indexing='ij')
                            block_csv = np.column_stack((lons.ravel(), lats.ravel(), block_val.astype(str).ravel()))
                            del lats, lons
                        del block, block_lat, block_lon
                        block_csv = block_csv[np.isfinite(block_val.ravel())]
                        del block_val
                        if block_csv.size > 0:
                            np.savetxt(f, block_csv, delimiter=delimiter, fmt='%s')
                        del block_csv
                        pbar.update(1)
        del grid

    def export_netcdf(self, data, name, caption='Exporting WGS84 NetCDF', engine='netcdf4', format='NETCDF3_64BIT'):
        """
        Export NetCDF file "velocity.nc":
        sbas.export_netcdf(velocity, 'velocity')
        """
        import xarray as xr
        import numpy as np
        import dask
        import os
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        from dask.distributed import utils_perf
        utils_perf.disable_gc_diagnosis()
    
        assert isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object'
    
        # convert the data to geographic coordinates if necessary
        if not self.is_geo(data):
            grid = self.ra2ll(data)
        else:
            grid = data
    
        filename = f'{name}.nc'
        if os.path.exists(filename):
            os.remove(filename)
        encoding = {data.name: self._compression(grid.shape)}
        delayed = grid.to_netcdf(filename, engine=engine, encoding=encoding, format=format, compute=False)
        tqdm_dask(dask.persist(delayed), desc=caption)
        del grid

    def export_vtk(self, data, name, caption='Exporting WGS84 VTK(s)', topo='auto', image=None, mask=None):
        """
        Export VTK file "velocity.vtk":
        sbas.export_vtk(velocity, 'velocity')
        """
        import xarray as xr
        import numpy as np
        from vtk import vtkStructuredGridWriter, vtkStringArray, VTK_BINARY
        from tqdm.auto import tqdm
        import os
    
        assert data is None or isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object or None'
        assert data is not None or image is not None, 'One of arguments "data" or "image" needs to be specified'
        
        # determine if data has a stack dimension and what it is
        stackvar = data.dims[0] if data is not None and len(data.dims) == 3 else None
        #print ('stackvar', stackvar)
        if stackvar is not None and np.issubdtype(data[stackvar].dtype, np.datetime64):
            stackvals = data[stackvar].dt.date.astype(str).values
        elif stackvar is not None:
            stackvals = data[stackvar].astype(str).values
        else:
            stackvals = [None]
        #print ('stackvals', stackvals)
    
        filename = f'{name}.vtk'
        if os.path.exists(filename):
            os.remove(filename)
    
        dss = []
        grid2d = None
    
        if data is not None:
            # convert the data to geographic coordinates if necessary
            data_ll = self.ra2ll(data) if not self.is_geo(data) else data
            # define 2D grid for interpolation 
            grid2d = data_ll.min(stackvar) if stackvar is not None else data_ll
            dss.append(data_ll)
        
        if image is not None:
            image_ll = self.ra2ll(image) if not self.is_geo(image) else image
            if not 'band' in image.dims:
                image_ll = image_ll.expand_dims('band')
            if grid2d is not None:
                image_ll = image_ll.interp_like(grid2d, method='linear')
                #.round().astype(np.uint8)
            else:
                grid2d = image_ll.isel(band=0)
            dss.append(image_ll.rename('colors'))
        
        #print ('dss', dss)
        #print ('grid2d', grid2d)
        
        if isinstance(topo, str) and topo == 'auto':
            dem = self.get_dem()
        elif topo is not None:
            # convert topography to geographic coordinates if necessary
            dem = self.ra2ll(topo) if not self.is_geo(topo) else topo
        if topo is not None:
            dem = dem.interp_like(grid2d, method='linear')
            if isinstance(mask, str) and mask == 'auto':
                dem = dem.where(np.isfinite(grid2d))
            elif mask is not None:
                dem = dem.where(mask.reindex_like(grid2d, method='nearest'))
            dss.append(dem.rename('z'))
    
        ds = xr.merge(dss, compat='override').rename({'lat': 'y', 'lon': 'x'})
        #print ('ds', ds)
    
        # prepare the progress bar
        with tqdm(desc=caption, total=len(stackvals)) as pbar:
            for stackidx, stackval in enumerate(stackvals):
                #print (stackidx, stackval)
                if stackval is not None:
                    out = ds.isel({stackvar: stackidx}).drop_vars(stackvar)
                    filename = f'{name}.{stackidx}.vtk'
                else:
                    out = ds
                    filename = f'{name}.vtk'
                #print ('out', out)
    
                vtk_grid = self.as_vtk(out)
                if stackval is not None:
                    metadata = vtkStringArray()
                    metadata.SetName(stackvar)
                    metadata.InsertNextValue(stackval)
                    vtk_grid.GetFieldData().AddArray(metadata)
    
                # convert to VTK structure and save to file
                writer = vtkStructuredGridWriter()
                writer.SetFileName(filename)
                writer.SetInputData(vtk_grid)
                writer.SetFileType(VTK_BINARY)
                writer.Write()
                pbar.update(1)
