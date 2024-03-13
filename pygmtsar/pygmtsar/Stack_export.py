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

    def export_stack_geotiff(self, data, name, caption='Exporting WGS84 GeoTIFFs'):
        """
        Export single GeoTIFF velocity.tif:
        sbas.export_stack(velocity_ps, 'velocity')
        
        Export pair-based stack of GeoTIFFs (like corr_ps.2024-01-08_2024-01-20.tif, ...):
        sbas.export_stack(corr_ps, 'corr_ps')
        
        Export date-based stack of GeoTIFFs (like disp_ps.2023-12-27.tif, ...):
        sbas.export_stack(disp_ps, 'disp_ps')

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

    def export_stack_geojson(self, data, name, caption='Exporting WGS84 GeoJSON', pivotal=True, digits=1):
        """
        Export single GeoJSON velocity.geojson:
        sbas.export_stack_geojson(velocity_ps, 'velocity')
    
        Export pair-based stack of GeoJSONs (like corr_ps.2024-01-08_2024-01-20.geojson, ...):
        sbas.export_stack_geojson(corr_ps, 'corr_ps')
    
        Export date-based stack of GeoJSONs (like disp_ps.2023-12-27.geojson, ...):
        sbas.export_stack_geojson(disp_ps, 'disp_ps')
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

