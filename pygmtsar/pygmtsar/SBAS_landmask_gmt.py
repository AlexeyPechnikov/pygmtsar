#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_merge import SBAS_merge

class SBAS_landmask_gmt(SBAS_merge):

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
