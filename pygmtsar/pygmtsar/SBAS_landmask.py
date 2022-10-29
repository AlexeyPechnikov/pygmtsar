#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_landmask_gmt import SBAS_landmask_gmt

class SBAS_landmask(SBAS_landmask_gmt):

    def set_landmask(self, landmask_filename):
        import os
        if landmask_filename is not None:
            self.landmask_filename = os.path.relpath(landmask_filename,'.')
        else:
            self.landmask_filename = None
        return self

    def get_landmask(self, inverse_geocode=False):
        import xarray as xr
        import os

        if self.landmask_filename is None:
            raise Exception('Set landmask first')

        # open DEM file and find the elevation variable
        # because sometimes grid includes 'crs' or other variables
        landmask = xr.open_dataset(self.landmask_filename, engine=self.engine, chunks=self.chunksize)
        assert 'lat' in landmask.coords and 'lon' in landmask.coords
        # define latlon array
        z_array_name = [data_var for data_var in landmask.data_vars if len(landmask.data_vars[data_var].coords)==2]
        assert len(z_array_name) == 1
        # extract the array and fill missed values by zero (mostly ocean area)
        landmask = landmask[z_array_name[0]].fillna(0)

        # select valid area only
        trans_dat = self.get_trans_dat()
        landmask = landmask.reindex_like(trans_dat, method='nearest')

        if inverse_geocode:
            return self.intf_ll2ra(landmask)
        return landmask
