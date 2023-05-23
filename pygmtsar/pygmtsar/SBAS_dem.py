#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_reframe import SBAS_reframe
from .PRM import PRM

class SBAS_dem(SBAS_reframe):

    def set_dem(self, dem_filename):
        import os
        if dem_filename is not None:
            self.dem_filename = os.path.relpath(dem_filename,'.')
        else:
            self.dem_filename = None
        return self

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small buffer produces incomplete area coverage and restricted NaNs
    # minimum buffer size: 8 arc seconds for 90 m DEM
    def get_dem(self, subswath=None, geoloc=False, buffer_degrees=0.02):
        import xarray as xr
        import os

        if self.dem_filename is None:
            raise Exception('Set DEM first')

        # open DEM file and find the elevation variable
        # because sometimes grid includes 'crs' or other variables
        dem = xr.open_dataset(self.dem_filename, engine=self.engine, chunks=self.chunksize)
        assert 'lat' in dem.coords and 'lon' in dem.coords, 'DEM should be defined as lat,lon grid'
        # define latlon array
        z_array_name = [data_var for data_var in dem.data_vars if len(dem.data_vars[data_var].coords)==2]
        assert len(z_array_name) == 1
        # extract the array and fill missed values (mostly water surfaces)
        dem = dem[z_array_name[0]].fillna(0)

        if geoloc is False:
            return dem

        bounds = self.get_master(subswath).dissolve().envelope.bounds.values[0].round(3)
        #print ('xmin, xmax', xmin, xmax)
        return dem\
                   .transpose('lat','lon')\
                   .sel(lat=slice(bounds[1]-buffer_degrees, bounds[3]+buffer_degrees),
                       lon=slice(bounds[0]-buffer_degrees, bounds[2]+buffer_degrees))

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small margin produces insufficient DEM not covers the defined area
    # https://docs.generic-mapping-tools.org/6.0/datasets/earth_relief.html
    # only bicubic interpolation supported as the best one for the case
    def download_dem(self, backend=None, product='SRTM1', resolution_meters=60, method=None, buffer_degrees=0.02, debug=False):
        """
        Use GMT server to download SRTM 1 or 3 arcsec data (@earth_relief_01s or @earth_relief_03s)
        Remove EGM96 geoid to make heights relative to WGS84
        Regrid to specified approximate resolution_meters (60m by default)
        """
        import pygmt
        import os
        #import subprocess
        from tqdm.auto import tqdm
        # 0.000833333333333 cell size for SRTM3 90m
        # 0.000277777777778 cell size for SRTM1 30m
        #scale = 0.000833333333333/90

        if self.dem_filename is not None:
            print ('NOTE: DEM exists, ignore the command. Use SBAS.set_dem(None) to allow new DEM downloading')
            return

        if backend is not None:
            print ('Note: backend argument is deprecated, just omit it')
        if method is not None:
            print ('Note: method argument is deprecated, just omit it')

        if product == 'SRTM1':
            resolution = '03s'
        elif product == 'SRTM3':
            resolution = '03s'
        elif product in ['01s', '03s']:
            pass
        else:
            print (f'ERROR: unknown product {product}. Available only SRTM1 ("01s") and SRTM3 ("03s") DEM using GMT servers')

        err, warn = self.validate()
        #print ('err, warn', err, warn)
        assert not err and not warn, 'ERROR: Please fix all the issues listed above to continue'

        # define approximate resolution in arc seconds
        spacing = f'{resolution_meters / 30}s'
        #print ('spacing', spacing)
        # generate DEM for the full area using GMT extent as W E S N
        minx, miny, maxx, maxy = self.df.dissolve().envelope.buffer(buffer_degrees).bounds.values[0]

        # Set the region for the grdcut and grdsample operations
        region = [minx, maxx, miny, maxy]

        gmtsar_sharedir = PRM().gmtsar_sharedir()
        geoid_filename = os.path.join(gmtsar_sharedir, 'geoid_egm96_icgem.grd')
        dem_filename = os.path.join(self.basedir, 'DEM_WGS84.nc')

        # use GMT commands pipeline to download and preprocess the DEM
        with tqdm(desc='DEM Downloading', total=1) as pbar:
            ortho = pygmt.datasets.load_earth_relief(resolution=resolution, region=region)
            ortho_resamp = pygmt.grdsample(ortho, region=region, spacing=spacing)
            geoid_resamp = pygmt.grdsample(geoid_filename, region=region, spacing=spacing)
            (ortho_resamp + geoid_resamp).to_netcdf(dem_filename)
    
            pbar.update(1)

        self.dem_filename = dem_filename
