#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_merge import SBAS_merge

class SBAS_landmask(SBAS_merge):

    def set_landmask(self, landmask_filename):
        import os
        if landmask_filename is not None:
            self.landmask_filename = os.path.relpath(landmask_filename,'.')
        else:
            self.landmask_filename = None
        return self

    def get_landmask(self, inverse_geocode=False, crop_valid=True):
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
        # round the coordinates up to 1 mm to have the same grid for dem and landmask
        landmask['lat'] = landmask.lat.round(8)
        landmask['lon'] = landmask.lon.round(8)

        if crop_valid:
            # select valid area only
            trans_dat = self.get_trans_dat()
            return landmask.reindex_like(trans_dat, method='nearest')

        if inverse_geocode:
            return self.intf_ll2ra(landmask)

        return landmask

    def download_landmask(self, backend=None, debug=False):
        """
        Use GMT local data or server to download and build landmask on interferogram DEM area
        """
        import pygmt
        import os
        from tqdm.auto import tqdm
        # correspond to SRTM3 DEM
        arcsec_degree = 0.000833333333333/3

        if self.landmask_filename is not None:
            print ('NOTE: landmask exists, ignore the command. Use SBAS.set_landmask(None) to allow new landmask downloading')
            return

        if backend is not None:
            print ('Note: backend argument is deprecated, just omit it')

        # generate the same as DEM grid
        landmask_filename = os.path.join(self.basedir, 'landmask.nc')

        dem = self.get_dem()
        scale = dem.lon.diff('lon')[0].item()
        llmin = dem.lon.min().item()
        llmax = dem.lon.max().item()
        ltmin = dem.lat.min().item()
        ltmax = dem.lat.max().item()
        region = f'{llmin}/{llmax}/{ltmin}/{ltmax}'
        #print('region', region)

        # define approximate resolution in arc seconds
        spacing = (dem.lat.diff('lat')[0]/arcsec_degree).round(3).item()
        # format as string
        spacing = f'{spacing}s'

        with tqdm(desc='Landmask Downloading', total=1) as pbar:
            landmask = pygmt.grdlandmask(resolution='f', region=region, spacing=spacing, maskvalues='NaN/1')
            if os.path.exists(landmask_filename):
                os.remove(landmask_filename)
            landmask.to_netcdf(landmask_filename)
            pbar.update(1)

        self.landmask_filename = landmask_filename
