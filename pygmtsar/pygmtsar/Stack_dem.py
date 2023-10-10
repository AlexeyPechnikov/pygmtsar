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

    # Buffer size in degrees to expand the area covered by the DEM
    buffer_degrees = 0.1
    
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
            self.dem_filename = os.path.relpath(dem_filename,'.')
        else:
            self.dem_filename = None
        return self

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small buffer produces incomplete area coverage and restricted NaNs
    # 0.02 degrees works well worldwide but not in Siberia
    # minimum buffer size: 8 arc seconds for 90 m DEM
    # subswath argument is required for aligning
    def get_dem(self, subswath=None, geoloc=False):
        """
        Retrieve the digital elevation model (DEM) data.

        Parameters
        ----------
        subswath : str, optional
            Subswath name. Default is None.
        geoloc : bool, optional
            Flag indicating whether to return geolocated DEM. If True, the returned DEM will be limited to the area covered
            by the specified subswath, plus an additional buffer around it. If False, the full DEM extent will be returned.
            Default is False.

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
        with zeros (mostly representing water surfaces). If the `geoloc` parameter is True, the returned DEM is geolocated,
        limited to the area covered by the specified subswath, plus an additional buffer around it. If `geoloc` is False,
        the full extent of the DEM will be returned.
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
        assert 'lat' in dem.coords and 'lon' in dem.coords, 'DEM should be defined as lat,lon grid'
        # define latlon array
        z_array_name = [data_var for data_var in dem.data_vars if len(dem.data_vars[data_var].coords)==2]
        assert len(z_array_name) == 1
        # extract the array and fill missed values (mostly water surfaces)
        dem = dem[z_array_name[0]].fillna(0)
        # round the coordinates up to 1 mm
        dem['lat'] = dem.lat.round(8)
        dem['lon'] = dem.lon.round(8)

        if geoloc is False:
            return dem

        bounds = self.get_reference(subswath).dissolve().envelope.bounds.values[0].round(3)
        #print ('xmin, xmax', xmin, xmax)
        return dem\
                   .transpose('lat','lon')\
                   .sel(lat=slice(bounds[1] - self.buffer_degrees, bounds[3] + self.buffer_degrees),
                        lon=slice(bounds[0] - self.buffer_degrees, bounds[2] + self.buffer_degrees))


    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # 0.02 degrees works well worldwide but not in Siberia
    # small margin produces insufficient DEM not covers the defined area
    # https://docs.generic-mapping-tools.org/6.0/datasets/earth_relief.html
    # only bicubic interpolation supported as the best one for the case
    def download_dem(self, product='1s', debug=False):
        """
        Download and preprocess digital elevation model (DEM) data.

        Parameters
        ----------
        product : str, optional
            Product type of the DEM data. Available options are '1s' or 'SRTM1' (1 arcsec ~= 30m, default)
            and '3s' or 'SRTM3' (3 arcsec ~= 90m).
        debug : bool, optional
            Enable debug mode. Default is False.

        Returns
        -------
        None

        Examples
        --------
        Download STRM1 DEM with a resolution of 30 meters and convert it to the default 60-meter grid:
        stack.download_dem()

        Download STRM1 DEM with a resolution of 30 meters and convert it to a 60-meter grid:
        stack.download_dem(resolution_meters=60)

        Download STRM3 DEM with a resolution of 90 meters and convert it to a 120-meter grid:
        stack.download_dem(product='STRM3', resolution_meters=120)

        Load and crop from local NetCDF file:
        stack.download_dem(product='GEBCO_2020/GEBCO_2020.nc')

        Load and crop from local GeoTIF file:
        stack.download_dem(product='GEBCO_2019.tif')

        Notes
        -----
        This method uses the GMT servers to download SRTM 1 or 3 arc-second DEM data. The downloaded data is then preprocessed
        by removing the EGM96 geoid to make the heights relative to the WGS84 ellipsoid. The DEM is regridded to the specified
        approximate resolution using bicubic interpolation.
        """
        import xarray as xr
        import numpy as np
        import pygmt
        # suppress warnings
        pygmt.config(GMT_VERBOSE='errors')
        import rioxarray as rio
        import os
        #import subprocess
        from tqdm.auto import tqdm

        if self.dem_filename is not None:
            print ('NOTE: DEM exists, ignore the command. Use Stack.set_dem(None) to allow new DEM downloading')
            return
        
        if product in ['SRTM1', '1s', '01s']:
            resolution = '01s'
        elif product in ['SRTM3', '3s', '03s']:
            resolution = '03s'
        else:
            # open from user-specified NetCDF or GeoTIF file
            resolution = None
            #print (f'ERROR: unknown product {product}. Available only SRTM1 ("01s") and SRTM3 ("03s") DEM using GMT servers')

        #err, warn = self.validate()
        #print ('err, warn', err, warn)
        #assert not err and not warn, 'ERROR: Please fix all the issues listed above to continue'

        # generate DEM for the area using GMT extent as W E S N
        # round the coordinates up to 1 mm
        minx, miny, maxx, maxy = self.df.dissolve().envelope.buffer(self.buffer_degrees).bounds.round(8).values[0]
        # Set the region for the grdcut and grdsample operations
        region = [minx, maxx, miny, maxy]

        gmtsar_sharedir = PRM().gmtsar_sharedir()
        geoid_filename = os.path.join(gmtsar_sharedir, 'geoid_egm96_icgem.grd')
        dem_filename = os.path.join(self.basedir, 'DEM_WGS84.nc')

        with tqdm(desc='DEM Downloading', total=1) as pbar:
            if resolution is not None:
                # download DEM
                ortho = pygmt.datasets.load_earth_relief(resolution=resolution, region=region)
            elif os.path.splitext(product)[-1] in ['.tiff', '.tif', '.TIF']:
                ortho = rio.open_rasterio(product, chunks=self.chunksize).squeeze(drop=True)\
                    .rename({'y': 'lat', 'x': 'lon'}).sel(lon=slice(minx, maxx))\
                    .drop('spatial_ref')
                if ortho.lat.diff('lat')[0].item() > 0:
                    ortho = ortho.sel(lat=slice(miny, maxy))
                else:
                    ortho = ortho.sel(lat=slice(maxy, miny))
                    ortho = ortho.reindex(lat=ortho.lat[::-1])
            elif os.path.splitext(product)[-1] in ['.nc', '.netcdf', '.grd']:
                ortho = xr.open_dataarray(product, engine=self.netcdf_engine, chunks=self.chunksize)\
                    .sel(lat=slice(miny, maxy), lon=slice(minx, maxx))

            #print ('ortho', ortho)
            geoid = xr.open_dataarray(geoid_filename).rename({'y': 'lat', 'x': 'lon'}).interp_like(ortho, method='cubic')
            #print ('geoid', geoid)
            if os.path.exists(dem_filename):
                os.remove(dem_filename)
            #print ('(ortho + geoid_resamp)', (ortho + geoid))
            encoding = {'dem': self.compression(ortho.shape)}
            (ortho + geoid).rename('dem').to_netcdf(dem_filename, encoding=encoding, engine=self.netcdf_engine)
            pbar.update(1)

        self.dem_filename = dem_filename
