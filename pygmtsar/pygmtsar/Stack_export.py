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
        Map a 2D image onto a 3D topography and convert it to a VTK structure for visualization.

        Parameters
        ----------
        dataset : xarray.Dataset
            The input dataset containing the 2D image and topography information.

        Returns
        -------
        vtk.vtkStructuredGrid
            A VTK structured grid representing the image mapped onto the topography.

        Notes
        -----
        This method converts an xarray dataset into a VTK structured grid, which can be used for 3D visualization. 
        The dataset should contain spatial dimensions and optionally topography information (e.g., elevation data).
        If the dataset contains RGB or RGBA bands, they will be properly handled and mapped as colors.
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

    def export_geotiff(self, data, name, caption='Exporting WGS84 GeoTIFF(s)', compress='LZW'):
        """
        Export the provided data to a GeoTIFF file.

        Parameters
        ----------
        data : xarray.DataArray
            The data to be exported as a GeoTIFF.
        name : str
            The base name for the GeoTIFF file(s).
        caption : str, optional
            A description for the export process, used for progress display. Default is 'Exporting WGS84 GeoTIFF(s)'.
        compress : str, optional
            The compression method to use for the GeoTIFF. Default is 'LZW'.

        Returns
        -------
        None
            The function writes the GeoTIFF file(s) to disk with the specified name.
            
        Examples
        --------
        Export a single GeoTIFF file "velocity.tif":
        sbas.export_geotiff(velocity, 'velocity')

        Export a stack of GeoTIFF files like corr.2024-01-01_2024-01-02.tif:
        sbas.export_geotiff(corr, 'corr')

        Export date-based stack of GeoTIFF files like disp.2024-01-01.tif:
        sbas.export_geotiff(disp, 'disp')
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
                self.as_geo(self.ra2ll(grid) if not self.is_geo(grid) else grid).rio.to_raster(filename, compress=compress)
                pbar.update(1)

    def export_geojson(self, data, name, caption='Exporting WGS84 GeoJSON', pivotal=True, digits=2, coord_digits=6):
        """
        Export the provided data to a GeoJSON file.

        Parameters
        ----------
        data : xarray.DataArray
            The data to be exported as GeoJSON.
        name : str
            The base name for the GeoJSON file(s).
        caption : str, optional
            A description for the export process, used for progress display. Default is 'Exporting WGS84 GeoJSON'.
        pivotal : bool, optional
            Whether to pivot the data. Default is True.
        digits : int, optional
            Number of decimal places to round the data values. Default is 2.
        coord_digits : int, optional
            Number of decimal places to round the coordinates. Default is 6.
    
        Returns
        -------
        None
            The function writes the GeoJSON file(s) to disk with the specified name.
            
        Examples
        --------
        Export a GeoJSON file "velocity.geojson":
        sbas.export_geojson(velocity, 'velocity')
        """
        import xarray as xr
        import geopandas as gpd
        import numpy as np
        from tqdm.auto import tqdm
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        try:
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()
        except ImportError:
            from distributed.gc import disable_gc_diagnosis
            disable_gc_diagnosis()

        assert isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object'

        # determine if data has a stack dimension and what it is
        stackvar = data.dims[0] if len(data.dims) == 3 else None

        # convert the data to geographic coordinates if necessary
        if not self.is_geo(data):
            grid = self.ra2ll(data)
        else:
            grid = data

        def block_as_json(block, stackvar, name):
            df = block.compute().to_dataframe().dropna().reset_index()
            df[name] = df[name].apply(lambda x: float(f"{x:.{digits}f}"))
            for col in df.columns:
                if np.issubdtype(df[col].dtype, np.datetime64):
                    df[col] = df[col].dt.date.astype(str)
            if stackvar is not None and np.issubdtype(df[stackvar].dtype, np.datetime64):
                df[stackvar] = df[stackvar].dt.date.astype(str)
            if stackvar is not None and pivotal:
                df = df.pivot_table(index=['lat', 'lon'], columns=stackvar,
                                    values=name, fill_value=np.nan).reset_index()
            # convert to geodataframe with value and geometry columns
            gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon.round(coord_digits), df.lat.round(coord_digits)))
            del df, gdf['lat'], gdf['lon']
            chunk_json = None
            if len(gdf):
                chunk_json = gdf.to_json(drop_id=True)
                # crop GeoJSON header and footer and split lines
                chunk_json = chunk_json[43:-2].replace('}}, ', '}},\n')
            del gdf
            return chunk_json

        # json header flag
        empty = True
        with open(f'{name}.geojson', 'w') as f:
            # GeoJSON header
            f.write('{"type": "FeatureCollection", "features": [\n')
    
            if 'stack' in data.dims:
                stack_blocks = np.array_split(np.arange(grid['stack'].size), np.arange(0, grid['stack'].size, self.chunksize1d)[1:])
                # prepare the progress bar
                with tqdm(desc=caption, total=len(stack_blocks)) as pbar:
                    for stack_block in stack_blocks:
                        block = grid.isel(stack=stack_block).drop_vars(['y','x'])
                        chunk_json = block_as_json(block, stackvar, data.name)
                        del block
                        if chunk_json is not None:
                            f.write(('' if empty else ',') + chunk_json)
                            empty = False
                        pbar.update(1)
            else:
                # split to equal chunks and rest
                # 1/4 NetCDF chunk is the smallest reasonable processing chunk
                lats_blocks = np.array_split(np.arange(grid.lat.size), np.arange(0, grid.lat.size, self.netcdf_chunksize//2)[1:])
                lons_blocks = np.array_split(np.arange(grid.lon.size), np.arange(0, grid.lon.size, self.netcdf_chunksize//2)[1:])
                # prepare the progress bar
                with tqdm(desc=caption, total=len(lats_blocks)*len(lons_blocks)) as pbar:
                    for lats_block in lats_blocks:
                        for lons_block in lons_blocks:
                            block = grid.isel(lat=lats_block, lon=lons_block)
                            chunk_json = block_as_json(block, stackvar, data.name)
                            del block
                            if chunk_json is not None:
                                f.write(('' if empty else ',') + chunk_json)
                                empty = False
                            pbar.update(1)

            # GeoJSON footer
            f.write(']}')
        del grid

    def export_csv(self, data, name, caption='Exporting WGS84 CSV', delimiter=',', digits=2, coord_digits=6):
        """
        Export the provided data to a CSV file.
    
        Parameters
        ----------
        data : xarray.DataArray
            The data to be exported as CSV.
        name : str
            The base name for the CSV file(s).
        caption : str, optional
            A description for the export process, used for progress display. Default is 'Exporting WGS84 CSV'.
        delimiter : str, optional
            The delimiter to use in the CSV file. Default is ','.
        digits : int, optional
            Number of decimal places to round the data values. Default is 2.
        coord_digits : int, optional
            Number of decimal places to round the coordinates. Default is 6.

        Returns
        -------
        None
            The function writes the CSV file(s) to disk with the specified name.
            
        Examples
        --------
        Export a CSV file "velocity.csv":
        sbas.export_csv(velocity, 'velocity')
        """
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        try:
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()
        except ImportError:
            from distributed.gc import disable_gc_diagnosis
            disable_gc_diagnosis()

        assert isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object'
    
        # convert the data to geographic coordinates if necessary
        if not self.is_geo(data):
            grid = self.ra2ll(data)
        else:
            grid = data

        # determine if data has a stack dimension and what it is
        if 'stack' in data.dims:
            stackvar = data.dims[0] if len(data.dims) == 2 else None
        else:
            stackvar = data.dims[0] if len(data.dims) == 3 else None
    
        with open(f'{name}.csv', 'w') as f:
            # CSV header
            f.write(delimiter.join(filter(None, [stackvar, 'lon', 'lat', data.name])) + '\n')
            if 'stack' in data.dims:
                stack_blocks = np.array_split(np.arange(grid['stack'].size), np.arange(0, grid['stack'].size, self.chunksize1d)[1:])
                # prepare the progress bar
                with tqdm(desc=caption, total=len(stack_blocks)) as pbar:
                    for stack_block in stack_blocks:
                        block = grid.isel(stack=stack_block).drop_vars(['y','x']).compute()
                        block_val = block.round(digits).values
                        block_lat = block.lat.round(coord_digits).values
                        block_lon = block.lon.round(coord_digits).values
                        if stackvar is not None:
                            stackvals = block[stackvar]
                            if np.issubdtype(stackvals.dtype, np.datetime64):
                                stackvals = stackvals.dt.date.astype(str)
                            block_csv = np.column_stack((np.repeat(stackvals, block_lon.size),
                                                         np.repeat(block_lon, stackvals.size),
                                                         np.repeat(block_lat, stackvals.size),
                                                         block_val.ravel()))
                        else:
                            block_csv = np.column_stack((block_lon, block_lat, block_val.astype(str).ravel()))
                        del block, block_lat, block_lon
                        block_csv = block_csv[np.isfinite(block_val.ravel())]
                        del block_val
                        if block_csv.size > 0:
                            np.savetxt(f, block_csv, delimiter=delimiter, fmt='%s')
                        del block_csv
                        pbar.update(1)
            else:
                # split to equal chunks and rest
                # 1/4 NetCDF chunk is the smallest reasonable processing chunk
                lats_blocks = np.array_split(np.arange(grid.lat.size), np.arange(0, grid.lat.size, self.netcdf_chunksize//2)[1:])
                lons_blocks = np.array_split(np.arange(grid.lon.size), np.arange(0, grid.lon.size, self.netcdf_chunksize//2)[1:])
                # prepare the progress bar
                with tqdm(desc=caption, total=len(lats_blocks)*len(lons_blocks)) as pbar:
                    for lats_block in lats_blocks:
                        for lons_block in lons_blocks:
                            block = grid.isel(lat=lats_block, lon=lons_block).compute()
                            block_val = block.round(digits).values
                            block_lat = block.lat.round(coord_digits).values
                            block_lon = block.lon.round(coord_digits).values
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
        Export the provided data to a NetCDF file.
    
        Parameters
        ----------
        data : xarray.DataArray
            The data to be exported as NetCDF.
        name : str
            The base name for the NetCDF file(s).
        caption : str, optional
            A description for the export process, used for progress display. Default is 'Exporting WGS84 NetCDF'.
        engine : str, optional
            The NetCDF engine to use (e.g., 'netcdf4'). Default is 'netcdf4'.
        format : str, optional
            The NetCDF format to use (e.g., 'NETCDF3_64BIT'). Default is 'NETCDF3_64BIT'.

        Returns
        -------
        None
            The function writes the NetCDF file to disk with the specified name.
            
        Examples
        --------
        Export a NetCDF file "velocity.nc":
        sbas.export_netcdf(velocity, 'velocity')
        """
        import xarray as xr
        import pandas as pd
        import numpy as np
        import dask
        import os
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        try:
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()
        except ImportError:
            from distributed.gc import disable_gc_diagnosis
            disable_gc_diagnosis()
    
        assert isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object'
    
        # convert the data to geographic coordinates if necessary
        if not self.is_geo(data):
            grid = self.ra2ll(data)
        else:
            grid = data

        if 'stack' in data.dims and isinstance(data.coords['stack'].to_index(), pd.MultiIndex):
            print (f"NOTE: open as xr.open_dataarray('{name}.nc').set_index(stack=['lat', 'lon'])")
            grid = grid.reset_index('stack')
    
        filename = f'{name}.nc'
        if os.path.exists(filename):
            os.remove(filename)
        encoding = {data.name: self._compression(grid.shape)}
        delayed = grid.to_netcdf(filename, engine=engine, encoding=encoding, format=format, compute=False)
        tqdm_dask(dask.persist(delayed), desc=caption)
        del grid

    def export_vtk(self, data, name, caption='Exporting WGS84 VTK(s)', topo='auto', image=None, mask=None):
        """
        Export the provided data to a VTK file.
    
        Parameters
        ----------
        data : xarray.DataArray or None
            The data to be exported as VTK. If None, `image` must be provided.
        name : str
            The base name for the VTK file(s).
        caption : str, optional
            A description for the export process, used for progress display. Default is 'Exporting WGS84 VTK(s)'.
        topo : str or None, optional
            The topography information to use for 3D mapping. If 'auto', the DEM is automatically retrieved. 
            If None, no topography will be applied. Default is 'auto'.
        image : xarray.DataArray or None, optional
            An optional image to overlay on the topography. Default is None.
        mask : xarray.DataArray or None, optional
            An optional mask to apply to the topography. If 'auto', the topography will be masked by the data. Default is None.
    
        Returns
        -------
        None
            The function writes the VTK file(s) to disk with the specified name.
    
        Examples
        --------
        Export a VTK file "velocity.vtk":
        sbas.export_vtk(velocity, 'velocity')
    
        Notes
        -----
        This function facilitates the visualization of InSAR data in VTK. By default, it includes topography data,
        and it can optionally overlay image data, such as satellite imagery, on the DEM to enhance the visualization.
        The `topo` option allows for automatic DEM inclusion, while the `mask` option helps manage NoData values,
        ensuring cleaner visual outputs.
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
