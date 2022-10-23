#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_reframe  import SBAS_reframe
from .PRM import PRM

class SBAS_dem_gmt_gdal(SBAS_reframe):

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small margin produces insufficient DEM not covers the defined area
    # https://docs.generic-mapping-tools.org/6.0/datasets/earth_relief.html
    def download_dem_gmt(self, product='SRTM1', resolution_meters=60, method='cubic', buffer_degrees=0.02, debug=False):
        """
        Use GMT server to download SRTM 1 or 3 arcsec data (@earth_relief_01s or @earth_relief_03s)
        Remove EGM96 geoid to make heights relative to WGS84
        Regrid to specified approximate resolution_meters (60m by default)
        """
        import os
        import subprocess
        from tqdm.auto import tqdm
        # 0.000833333333333 cell size for SRTM3 90m
        # 0.000277777777778 cell size for SRTM1 30m
        scale = 0.000833333333333/90

        if self.dem_filename is not None:
            print ('NOTE: DEM exists, ignore the command. Use SBAS.set_dem(None) to allow new DEM downloading')
            return

        if product == 'SRTM1':
            gridname = '@earth_relief_01s'
        elif product == 'SRTM3':
            gridname = '@earth_relief_03s'
        else:
            print (f'ERROR: unknown product {product}. Available only SRTM1 and SRTM3 DEM using GMT servers')
    
        if method != 'cubic':
            print ('NOTE: only bicubic interpolation supported as the best one for the case')
    
        err, warn = self.validate()
        #print ('err, warn', err, warn)
        assert not err and not warn, 'ERROR: Please fix all the issues listed above to continue'

        # define resolution in decimal degrees
        resolution_degrees = scale * resolution_meters
        # generate DEM for the full area using GMT extent as W E S N
        minx, miny, maxx, maxy = self.df.dissolve().envelope.buffer(buffer_degrees).bounds.values[0]
    
        gmtsar_sharedir = PRM().gmtsar_sharedir()
        geoid_filename = os.path.join(gmtsar_sharedir, 'geoid_egm96_icgem.grd')
        ortho_filename = os.path.join(self.basedir, 'dem_ortho.grd')
        ortho_resamp_filename = os.path.join(self.basedir, 'dem_ortho_resamp.grd')
        geoid_resamp_filename = os.path.join(self.basedir, 'geoid_resamp.grd')
        dem_filename = os.path.join(self.basedir, 'DEM_WGS84.nc')

        # helper function to run external commands
        def run_cmd(argv):
            if debug:
                print ('DEBUG: argv', argv)
            cmd = subprocess.run(argv, capture_output=True)
            if cmd.returncode != 0:
                print (cmd.stderr)
                raise Exception('DEM processing error using GMT backend')
            if debug:
                print (cmd.stdout)

        # use GMT commands pipeline to download and preprocess the DEM
        with tqdm(desc='DEM Downloading', total=1) as pbar:
            # get srtm data
            argv = ['gmt', 'grdcut', gridname, f'-R{minx}/{maxx}/{miny}/{maxy}', f'-G{ortho_filename}', '-V']
            run_cmd(argv)
            argv = ['gmt', 'grdsample', f'-I{resolution_degrees}', ortho_filename, f'-G{ortho_resamp_filename}', '-V']
            run_cmd(argv)
            # resample and remove geoid
            argv = ['gmt', 'grdsample', geoid_filename, f'-R{ortho_resamp_filename}', f'-G{geoid_resamp_filename}', '-V']
            run_cmd(argv)
            argv = ['gmt', 'grdmath', '-V', ortho_resamp_filename, geoid_resamp_filename, 'ADD', '=', dem_filename]
            run_cmd(argv)
            pbar.update(1)
    
        # cleanup
        for filename in [ortho_resamp_filename, ortho_filename, geoid_resamp_filename]:
            if os.path.exists(filename):
                os.remove(filename)
    
        self.dem_filename = dem_filename

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small margin produces insufficient DEM not covers the defined area
    def download_dem_old(self, product='SRTM1', resolution_meters=60, method='cubic', buffer_degrees=0.02, debug=False):
        """
        Download SRTM3 90m or SRTM1 30 m DEM
        Regrid the DEM to specified approximate resolution_meters 60m by default
        Use resampling method 'cubic' by default (see gdalwarp documentation for -r <resampling_method>)
        """
        import urllib.request
        import elevation
        import os
        import subprocess
        from tqdm.auto import tqdm
        import joblib
        # 0.000833333333333 cell size for SRTM3 90m
        # 0.000277777777778 cell size for SRTM1 30m
        scale = 0.000833333333333/90
        
        if self.dem_filename is not None:
            print ('NOTE: DEM exists, ignore the command. Use SBAS.set_dem(None) to allow new DEM downloading')
            return
        
        resolution_degrees = scale * resolution_meters

        err, warn = self.validate()
        #print ('err, warn', err, warn)
        assert not err and not warn, 'ERROR: Please fix all the issues listed above to continue'

        gtx_url = 'https://github.com/mobigroup/proj-datumgrid/blob/master/egm96_15.gtx?raw=true'
        gtx_filename = os.path.join(self.basedir, 'egm96_15.gtx')
        tif_filename = os.path.join(self.basedir, 'DEM_EGM96.tif')
        grd_filename = os.path.join(self.basedir, 'DEM_WGS84.nc')

        if not os.path.exists(gtx_filename):
            with urllib.request.urlopen(gtx_url) as fin:
                with open(gtx_filename, 'wb') as fout:
                    fout.write(fin.read())

        #if not os.path.exists(tif_filename):
        # generate DEM for the full area
        bounds = self.df.dissolve().envelope.bounds.values[0]
        # show progress indicator
        def func():
            # left bottom right top
            elevation.clip(bounds=bounds,
                            product=product,
                            margin=str(buffer_degrees),
                            output=os.path.realpath(tif_filename))
            elevation.clean()
        with self.tqdm_joblib(tqdm(desc='Downloading', total=1)) as progress_bar:
            _ = joblib.Parallel(n_jobs=1)(joblib.delayed(func)() for i in [1])

        #if not os.path.exists(grd_filename):
        # convert to WGS84 ellipsoidal heights
        argv = ['gdalwarp', '-co', 'COMPRESS=DEFLATE',
                '-r', method,
                '-s_srs', f'+proj=longlat +datum=WGS84 +no_defs +geoidgrids=./egm96_15.gtx',
                '-t_srs', '+proj=longlat +datum=WGS84 +no_def',
                '-tr', str(resolution_degrees), str(resolution_degrees),
                '-overwrite',
                '-ot', 'Float32', '-of', 'NetCDF',
                'DEM_EGM96.tif', 'DEM_WGS84.nc']
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: download_dem', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: download_dem', stdout_data.decode('ascii'))

        self.dem_filename = grd_filename
