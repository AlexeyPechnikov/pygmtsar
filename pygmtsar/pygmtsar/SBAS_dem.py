#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_dem_gmt_gdal  import SBAS_dem_gmt_gdal

class SBAS_dem(SBAS_dem_gmt_gdal):

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

    # wrapper
    def download_dem(self, backend=None, **kwargs):
        if backend is None:
            return self.download_dem_old(**kwargs)
        elif backend == 'GMT':
            return self.download_dem_gmt(**kwargs)
        else:
            raise Exception(f'Unknown backend {backend}. Use None or GMT')
