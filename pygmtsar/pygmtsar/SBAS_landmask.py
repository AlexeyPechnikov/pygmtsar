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

    def get_landmask(self, subswath=None, geoloc=False, buffer_degrees=0.02, inverse_geocode=False):
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

        if inverse_geocode:
            return self.intf_ll2ra(subswath, landmask)

        dem = self.get_dem(subswath, geoloc, buffer_degrees)
        return landmask.reindex_like(dem, method='nearest')

    def download_landmask(self, backend='GMT', debug=False):
        """
        Use GMT local data or server to download and build landmask on interferogram DEM area
        """
        import os
        import subprocess
        from tqdm.auto import tqdm

        if self.landmask_filename is not None:
            print ('NOTE: landmask exists, ignore the command. Use SBAS.set_landmask(None) to allow new landmask downloading')
            return

        # generate the same as DEM grid
        landmask_filename = os.path.join(self.basedir, 'landmask.nc')

        dem = self.get_dem()
        scale = dem.lon.diff('lon')[0].item()
        llmin = dem.lon.min().item()
        llmax = dem.lon.max().item()
        ltmin = dem.lat.min().item()
        ltmax = dem.lat.max().item()

        # helper function to run external commands
        def run_cmd(argv):
            if debug:
                print ('DEBUG: argv', argv)
            cmd = subprocess.run(argv, capture_output=True)
            if cmd.returncode != 0:
                print (cmd.stderr)
                raise Exception('DEM processing error using GMT backend')
            if debug:
                print (cmd.stderr)
                print (cmd.stdout)

        # use GMT commands pipeline to download and preprocess the DEM
        with tqdm(desc='Landmask Downloading', total=1) as pbar:
            argv = ['gmt', 'grdlandmask', f'-G{landmask_filename}',
                    f'-R{llmin}/{llmax}/{ltmin}/{ltmax}', f'-I{scale}', '-V', '-NNaN/1', '-Df']
            run_cmd(argv)
            pbar.update(1)

        self.landmask_filename = landmask_filename
