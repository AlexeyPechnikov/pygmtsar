# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_reframe import Stack_reframe
from .PRM import PRM

class Stack_dem(Stack_reframe):

    buffer_degrees = 0.02

    def get_extent_ra(self):
        """
        minx, miny, maxx, maxy = np.round(geom.bounds).astype(int)
        """
        import numpy as np
        from shapely.geometry import LineString

        dem = self.get_dem()
        df = dem.isel(lon=[0,-1]).to_dataframe().reset_index()
        geom = self.geocode(LineString(np.column_stack([df.lon, df.lat])))
        return geom

#     def get_extent(self, grid=None, subswath=None):
#         import numpy as np
# 
#         bounds = self.get_bounds(self.get_reference(subswath))
#         if grid is None:
#             return bounds
#         return grid\
#                .transpose('lat','lon')\
#                .sel(lat=slice(bounds[1], bounds[3]),
#                     lon=slice(bounds[0], bounds[2]))

    def get_geoid(self, grid=None):
        """
        Get EGM96 geoid heights.

        Parameters
        ----------
        grid : xarray array or dataset, optional
            Interpolate geoid heights on the grid. Default is None.

        Returns
        -------
        None

        Examples
        --------
        stack.get_geoid()

        Notes
        -----
        See EGM96 geoid heights on http://icgem.gfz-potsdam.de/tom_longtime
        """
        import xarray as xr
        import os

        gmtsar_sharedir = PRM().gmtsar_sharedir()
        geoid_filename = os.path.join(gmtsar_sharedir, 'geoid_egm96_icgem.grd')
        geoid = xr.open_dataarray(geoid_filename, engine=self.netcdf_engine, chunks=self.netcdf_chunksize).rename({'y': 'lat', 'x': 'lon'})
        if grid is not None:
            geoid = geoid.interp_like(grid, method='cubic')
        return geoid

    def set_dem(self, dem_filename):
        """
        Set the filename of the digital elevation model (DEM) in WGS84 NetCDF grid.

        Parameters
        ----------
        dem_filename : str or None
            Path to the WGS84 NetCDF DEM file. If provided, the DEM filename will be set. If None, the DEM filename will be cleared.

        Returns
        -------
        self
            Returns the modified instance of the class.

        Examples
        --------
        Set the DEM filename:
        stack = stack.set_dem('data/DEM_WGS84.nc')

        Notes
        -----
        This method sets the filename of the digital elevation model (DEM) data to be used in the SAR processing.
        The DEM file is a NetCDF dataset in geographical coordinates containing elevation values for the Earth's surface.
        By setting the DEM filename, it allows SAR processing algorithms to utilize the DEM data for geocoding, topographic
        correction, and other applications. If `dem_filename` is None, the DEM filename will be cleared.
        """
        import os
        if dem_filename is not None:
            assert os.path.exists(dem_filename), f'DEM file not found: {dem_filename}'
            assert os.path.isfile(dem_filename) and os.access(dem_filename, os.R_OK), f'DEM file is not readable: {dem_filename}'
            self.dem_filename = os.path.relpath(dem_filename,'.')
        else:
            self.dem_filename = None
        return self

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small buffer produces incomplete area coverage and restricted NaNs
    # 0.02 degrees works well worldwide but not in Siberia
    # minimum buffer size: 8 arc seconds for 90 m DEM
    # subswath argument is required for aligning
    def get_dem(self, subswath=None):
        """
        Retrieve the digital elevation model (DEM) data.

        Parameters
        ----------
        subswath : str, optional
            Subswath name. Default is None.

        Returns
        -------
        xarray.DataArray
            The DEM data as a DataArray object.

        Raises
        ------
        Exception
            If the DEM is not set before calling this method.

        Examples
        --------
        Get DEM for all the processed subswaths:
        topo_ll = stack.get_dem()

        Get DEM for a single subswath IW1:
        topo_ll = stack.get_dem(1)

        Notes
        -----
        This method retrieves the digital elevation model (DEM) data previously downloaded and stored in a NetCDF file.
        The DEM file is opened, and the elevation variable is extracted. Any missing values in the elevation data are filled
        with zeros (mostly representing water surfaces).
        """
        import xarray as xr
        import os
        import warnings
        # supress warnings "UserWarning: The specified chunks separate the stored chunks along dimension"
        warnings.filterwarnings('ignore')

        if self.dem_filename is None:
            raise Exception('Set DEM first')

        # open DEM file and find the elevation variable
        # because sometimes grid includes 'crs' or other variables
        dem = xr.open_dataset(self.dem_filename, engine=self.netcdf_engine, chunks=self.chunksize)
        if 'lat' not in dem.coords and 'y' in dem.coords:
            dem = dem.rename({'y': 'lat', 'x': 'lon'})
        # define latlon array
        z_array_name = [data_var for data_var in dem.data_vars if len(dem.data_vars[data_var].coords)==2]
        assert len(z_array_name) == 1
        # extract the array and fill missed values (mostly water surfaces)
        dem = dem[z_array_name[0]]
        #.fillna(0)
        # round the coordinates up to 1 mm
        dem['lat'] = dem.lat.round(8)
        dem['lon'] = dem.lon.round(8)

        # crop to reference scene and subswath if defined
        bounds = self.get_bounds(self.get_reference(subswath))
        return dem\
               .transpose('lat','lon')\
               .sel(lat=slice(bounds[1] - self.buffer_degrees, bounds[3] + self.buffer_degrees),
                    lon=slice(bounds[0] - self.buffer_degrees, bounds[2] + self.buffer_degrees))
        #return self.get_extent(dem, subswath)

    def load_dem(self, data, geometry='auto', buffer_degrees=None):
        """
        Load and preprocess digital elevation model (DEM) data from specified datafile or variable.

        Parameters
        ----------
        data : xarray dataarray or str
            DEM filename or variable.

        Returns
        -------
        None

        Examples
        --------
        Load and crop from local NetCDF file:
        stack.load_dem('GEBCO_2020/GEBCO_2020.nc')

        Load and crop from local GeoTIF file:
        stack.load_dem('GEBCO_2019.tif')

        Load from Xarray DataArray or Dataset:
        stack.set_dem(None).load_dem(dem)
        stack.set_dem(None).load_dem(dem.to_dataset())

        Notes
        -----
        This method loads DEM from the user specified file. The data is then preprocessed by removing
        the EGM96 geoid to make the heights relative to the WGS84 ellipsoid.
        """
        import xarray as xr
        import numpy as np
        import rioxarray as rio
        import geopandas as gpd
        import os

        dem_filename = os.path.join(self.basedir, 'DEM_WGS84.nc')

        if self.dem_filename is not None:
            print ('NOTE: DEM exists, ignore the command. Use Stack.set_dem(None) to allow new DEM downloading')
            return

        if buffer_degrees is not None:
            self.buffer_degrees = buffer_degrees

        if isinstance(data, (xr.Dataset)):
            ortho = data[list(data.data_vars)[0]]
        elif isinstance(data, (xr.DataArray)):
            ortho = data
        elif isinstance(data, str) and os.path.splitext(data)[-1] in ['.tiff', '.tif', '.TIF']:
            ortho = rio.open_rasterio(data, chunks=self.chunksize).squeeze(drop=True)\
                .rename({'y': 'lat', 'x': 'lon'})\
                .drop('spatial_ref')
            if ortho.lat.diff('lat')[0].item() < 0:
                ortho = ortho.reindex(lat=ortho.lat[::-1])
        elif isinstance(data, str) and os.path.splitext(data)[-1] in ['.nc', '.netcdf', '.grd']:
            ortho = xr.open_dataarray(data, engine=self.netcdf_engine, chunks=self.chunksize)
        elif isinstance(data, str):
            print ('ERROR: filename extension is not recognized. Should be one from .tiff, .tif, .TIF, .nc, .netcdf, .grd')
        else:
            print ('ERROR: argument is not an Xarray object and it is not a file name')

        # crop to reference scene extent
        geometry_auto = type(geometry) == str and geometry == 'auto'
        bounds = self.get_bounds(self.get_reference() if geometry_auto else geometry)
        ortho = ortho\
               .transpose('lat','lon')\
               .sel(lat=slice(bounds[1] - self.buffer_degrees, bounds[3] + self.buffer_degrees),
                    lon=slice(bounds[0] - self.buffer_degrees, bounds[2] + self.buffer_degrees))

        # heights correction
        geoid = self.get_geoid(ortho)
        if os.path.exists(dem_filename):
            os.remove(dem_filename)
        encoding = {'dem': self._compression(ortho.shape)}
        (ortho + geoid).rename('dem').load()\
            .to_netcdf(dem_filename, encoding=encoding, engine=self.netcdf_engine)

        self.dem_filename = dem_filename

    def download_dem(self, geometry='auto', product='1s'):
        print ('NOTE: Function is removed. Download DEM using Tiles().download_dem()')
        print ('and load with Stack.load_dem() function.')
