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
from .NCubeVTK import NCubeVTK

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

    def export_vtk(self, data, name, caption='Exporting WGS84 VTK(s)', topo='auto', image=None, band_mask=None, use_sealevel=False):
        """
        Export VTK file "velocity.vtk":
        sbas.export_vtk(velocity, 'velocity')
        """
        import xarray as xr
        import numpy as np
        from vtk import vtkStructuredGridWriter, vtkStringArray
        from tqdm.auto import tqdm
        import os
    
        assert isinstance(data, xr.DataArray), 'Argument data is not an xr.DataArray object'
    
        # determine if data has a stack dimension and what it is
        stackvar = data.dims[0] if len(data.dims) == 3 else None
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
    
        # convert the data to geographic coordinates if necessary
        if not self.is_geo(data):
            data_ll = self.ra2ll(data)
        else:
            data_ll = data
        # define 2D grid for interpolation 
        grid2d = data_ll.min(stackvar) if stackvar is not None else data_ll
        #print ('grid2d', grid2d)
        
        dss = []
        if isinstance(topo, str) and topo == 'auto':
            dem = self.get_dem()
        elif topo is not None:
            # convert topography to geographic coordinates if necessary
            dem = self.ra2ll(topo) if not self.is_geo(topo) else topo
        if topo is not None:
            dem = dem.reindex_like(grid2d, method='nearest').where(np.isfinite(grid2d))
            dss.append(dem.rename('z'))
    
        if image is not None:
            dss.append(image.reindex_like(grid2d, method='nearest').round().astype(np.uint8).rename('colors'))
        #print ('dss', dss)
    
        # prepare the progress bar
        with tqdm(desc=caption, total=len(stackvals)) as pbar:
            for stackidx, stackval in enumerate(stackvals):
                #print (stackidx, stackval)
                if stackval is not None:
                    ds = xr.merge([*dss, data_ll.isel({stackvar: stackidx})], compat='override')\
                        .rename({'lat': 'y', 'lon': 'x'})\
                        .drop_vars(stackvar)
                    filename = f'{name}.{stackidx}.vtk'
                else:
                    ds = xr.merge([*dss, data_ll], compat='override').rename({'lat': 'y', 'lon': 'x'})
                    filename = f'{name}.vtk'
                #print ('ds', ds)
                
                vtk_grid = NCubeVTK.ImageOnTopography(ds)
                if stackval is not None:
                    metadata = vtkStringArray()
                    metadata.SetName(stackvar)
                    metadata.InsertNextValue(stackval)
                    vtk_grid.GetFieldData().AddArray(metadata)

                # convert to VTK structure and save to file
                writer = vtkStructuredGridWriter()
                writer.SetFileName(filename)
                writer.SetInputData(vtk_grid)
                writer.Write()
                pbar.update(1)
