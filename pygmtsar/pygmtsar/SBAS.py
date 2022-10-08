#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to scan *xml files and orbits and make data.in file like to prep_data_linux.csh & prep_data.csh tools
#import pytest

class SBAS:

    # NetCDF options
    netcdf_chunksize = 512
    netcdf_compression = dict(zlib=True, complevel=3, chunksizes=(netcdf_chunksize,netcdf_chunksize))
    netcdf_engine = 'h5netcdf'
    
    import contextlib
    @staticmethod
    @contextlib.contextmanager
    def tqdm_joblib(tqdm_object):
        import joblib
        """Context manager to patch joblib to report into tqdm progress bar given as argument"""
        class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)

            def __call__(self, *args, **kwargs):
                tqdm_object.update(n=self.batch_size)
                return super().__call__(*args, **kwargs)

        old_batch_callback = joblib.parallel.BatchCompletionCallBack
        joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
        try:
            yield tqdm_object
        finally:
            joblib.parallel.BatchCompletionCallBack = old_batch_callback
            tqdm_object.close()

    #text2date('V20171110T225942'), text2date('20171117t145927')
    @staticmethod
    def text2date(text, as_date=True):
        from datetime import datetime
        
        date_fmt = '%Y%m%dT%H%M%S'
        date_str = text.replace('V','').replace('t','T')
        dt = datetime.strptime(date_str, date_fmt)
        if as_date == False:
            return dt
        return dt.date()

    @staticmethod
    def is_ra(grid):
        dims = grid.dims
        if 'y' in dims and 'x' in dims:
            return True
        return False

    @staticmethod
    def is_geo(grid):
        dims = grid.dims
        if 'lat' in dims and 'lon' in dims:
            return True
        return False

    def as_geo(self, da):
        """
        Add spatial attributes to allow use rioxarray functions like to .rio.clip([geometry])
        """
        import sys
        assert 'rioxarray' in sys.modules, 'rioxarray module is not found'
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

    @staticmethod
    def is_same(grid1, grid2):
        dims1 = grid1.dims
        dims2 = grid2.dims
        if 'lat' in dims1 and 'lon' in dims1 and 'lat' in dims2 and 'lon' in dims2:
            return True
        if 'y' in dims1 and 'x' in dims1 and 'y' in dims2 and 'x' in dims2:
            return True
        return False

    @staticmethod
    def nearest_grid(in_grid, search_radius_pixels=300):
        from pygmtsar import PRM
        return PRM.nearest_grid(in_grid, search_radius_pixels)

    def snaphu_config(self, defomax=0, **kwargs):
        return self.PRM().snaphu_config(defomax, **kwargs)
 
    def cropna(self, das):
        # crop NaNs
        dims = [dim for dim in das.dims if dim != 'pair' and dim != 'date']
        dim0 = [dim for dim in das.dims if dim in ['pair', 'date']]
        assert len(dims) == 2, 'ERROR: interferogram should be 2D array'
        da = das.min(dim0)
        indexer = {}
        for dim in dims:
            da = da.dropna(dim=dim, how='all')
            dim_min, dim_max = da[dim].min().item(), da[dim].max().item()
            indexer[dim] = slice(dim_min, dim_max)
        #print ('indexer', indexer)
        return das.loc[indexer]
    
    def gaussian(self, da, wavelength, truncate=3.0, approximate=False, debug=False):
        """
        Gaussian filter for arrays with NaN values.
            da - input dataarray with NaNs allowed,
            wavelength - cut-off wavelength [m],
            truncate - filter window size [sigma],
            approximate - boolean flag to use approximate calculation based on filtering on decimated array,
            debug - pront debug information.
        Returns filtered dataarray with the same coordinates as input one.
        Fast approximate calculation silently skipped when sigma is less than 64 so the result is always exact for small filters.
        """
        import xarray as xr
        import numpy as np
        from pygmtsar import PRM

        # ground pixel size
        dy, dx = self.pixel_size(da)
        # gaussian kernel
        sigma_y = np.round(wavelength / dy)
        sigma_x = np.round(wavelength / dx)
    
        if debug:
            print ('DEBUG: Gaussian filtering sigma_y, sigma_x', sigma_y, sigma_x)
    
        # there is a compromise between accuracy and speed
        dec = 32
        # ignore decimation when it is less than or equal to 1
        dec_y = int(np.floor(sigma_y/dec))
        dec_x = int(np.floor(sigma_x/dec))
        
        if approximate and dec_y > 1 and dec_x > 1:
            # filter decimated array for faster and less accurate processing
            if debug:
                print ('DEBUG: Gaussian filtering decimation dec_y, dec_x', dec_y, dec_x)
        
            # decimate and filter
            da_dec = da.coarsen({'y': dec_y, 'x': dec_x}, boundary='pad').mean()
            values  = PRM.nanconvolve2d_gaussian(da_dec.values,
                                                      (sigma_y/dec_y, sigma_x/dec_x),
                                                      truncate=truncate)
            da_dec = xr.DataArray(values, coords=da_dec.coords)

            # scale to original grid and filter with small kerkel equal to the decimation
            da_dec = da_dec.interp_like(da, method='nearest')
            values   = PRM.nanconvolve2d_gaussian(da_dec.values,
                                                  (dec_y, dec_x),
                                                  truncate=truncate)
        else:
            values   = PRM.nanconvolve2d_gaussian(da.values, (sigma_y,sigma_x), truncate=truncate)
        return xr.DataArray(values, coords=da.coords)

    def __repr__(self):
        return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    def dump(self, to_path=None):
        import pickle
        import os

        if to_path is None:
            sbas_pickle = os.path.join(self.basedir, 'SBAS.pickle')
        else:
            if os.path.isdir(to_path):
                sbas_pickle = os.path.join(to_path, 'SBAS.pickle')
            else:
                sbas_pickle = to_path
    
        print (f'NOTE: save state to file {sbas_pickle}')
        pickle.dump(self, open(sbas_pickle, 'wb'))

        return

    @staticmethod
    def restore(from_path):
        import pickle
        import os

        if os.path.isdir(from_path):
            sbas_pickle = os.path.join(from_path, 'SBAS.pickle')
        else:
            sbas_pickle = from_path

        print (f'NOTE: load state from file {sbas_pickle}')
        return pickle.load(open(sbas_pickle, 'rb'))

    def __init__(self, datadir, dem_filename=None, basedir='.', landmask_filename=None,
                filter_orbit=None, filter_mission=None, filter_subswath=None, filter_polarization=None, force=True):
        import os
        import shutil
        from glob import glob
        import pandas as pd
        import geopandas as gpd
        import shapely
        import numpy as np
        from datetime import datetime
        from dateutil.relativedelta import relativedelta
        oneday = relativedelta(days=1)

        def pattern2paths(pattern):
            path_pattern = os.path.join(datadir, '**', pattern)
            paths = glob(path_pattern, recursive=True)
            return paths

        assert filter_orbit is None or filter_orbit=='A' or filter_orbit=='D', \
            'ERROR: use symbol A (Ascending) or D (Descending) for orbit filter'
        assert filter_mission is None or filter_mission=='S1A' or filter_mission=='S1B', \
            'ERROR: use name S1A or S1B for mission filter'
        assert filter_subswath is None or filter_subswath in [1,2,3,12,23,123], \
            'ERROR: use a single or sequential numbers 1, 2, 3, 12, 23, 123 for subswath filter'
        assert filter_polarization is None or filter_polarization in ['VV','VH','HH','HV'], \
            'ERROR: use VV or VH or HH or HV for polarization filter'

        # processing directory
        if basedir is None:
            self.basedir = '.'
        else:
            # (re)create basedir only when force=True
            if force and os.path.exists(basedir):
                shutil.rmtree(basedir)
            os.makedirs(basedir, exist_ok=True)
            self.basedir = basedir

        if filter_polarization is None:
            filter_polarization = '??'
        if filter_subswath is None:
            filter_subswath  = '?'
        else:
            filter_subswath = f'[{filter_subswath}]'
        # filter mission
        if filter_mission is not None:
            path_pattern = f'{filter_mission.lower()}-iw{filter_subswath}-slc-{filter_polarization.lower()}-*'
        else:
            path_pattern = f's1?-iw{filter_subswath}-slc-{filter_polarization.lower()}-*'
        datapaths = pattern2paths(path_pattern + '.tiff')
        #print ('datapaths', datapaths)
        metapaths = pattern2paths(path_pattern + '.xml')
        #print ('metapaths', metapaths)

        datanames = [os.path.splitext(os.path.split(path)[-1])[0] for path in datapaths]
        #print ('datanames', datanames)
        metanames = [os.path.splitext(os.path.split(path)[-1])[0] for path in metapaths]
        #print ('metanames', metanames)

        datas = dict(zip(datanames, datapaths))
        metas = dict(zip(metanames, metapaths))

        # define the same order when and only when the names are the same
        datanames = sorted(datanames)
        metanames = sorted(metanames)
        assert datanames == metanames, 'Found inconsistent set of .tiff and .xml files'
        # reorder paths using the same order
        datapaths = [datas[name] for name in datanames]
        metapaths = [metas[name] for name in metanames]

        # points to datadir and extensions tiff, xml
        #print ('filenames', filenames)
        dts = [self.text2date(name.split('-')[4],False) for name in datanames]
        #print ('filedatetimes', dts)

        ds = [dt.date() for dt in dts]
        #print ('filedates', ds)

        df = pd.DataFrame({'date':[str(d) for d in ds], 'datetime': dts, 'datapath': datapaths, 'metapath': metapaths})
        #print ('self.df', self.df)

        # filter mission and always ignore approximate RESORB orbits to download precise POEORB when possible
        if filter_mission is not None:
            orbit_path_pattern = f'{filter_mission.upper()}_OPER_AUX_POEORB_OPOD_*.EOF'
        else:
            orbit_path_pattern = 'S1?_OPER_AUX_POEORB_OPOD_*.EOF'
        orbitpaths = pattern2paths(orbit_path_pattern)
        #print ('orbitpaths', orbitpaths)
        orbitnames = [os.path.splitext(os.path.split(path)[-1])[0] for path in orbitpaths]
        if orbitpaths:
            orbit_dates = [(self.text2date(name.split('_')[-2]), self.text2date(name.split('_')[-1])) for name in orbitnames]
            orbits = dict(zip(orbit_dates, orbitpaths))
            #print ('orbits', orbits)
            # look for as precise (from date-1 day to date+1 day) as restituted orbits (from date to date)
            orbits = [orbits.get((date-oneday, date+oneday)) or orbits.get((date,date)) for date in ds]
            #print ('fileorbits', fileorbits)
            df['orbitpath'] = orbits
        else:
            df['orbitpath'] = None

        # add some calculated properties
        df['subswath'] = [int(filename.split('-')[1][-1]) for filename in datanames]
        df['mission'] = [filename.split('-')[0].upper() for filename in datanames]
        df['polarization'] = [filename.split('-')[3].upper() for filename in datanames]
    
        # read approximate locations
        geolocs = [shapely.geometry.MultiPoint(self.geoloc(path).geometry).minimum_rotated_rectangle for path in metapaths]
        #print ('geolocs', geolocs)
        df = gpd.GeoDataFrame(df, geometry=geolocs)

        # define orbit directions
        orbits = [self.annotation(path)['product']['generalAnnotation']['productInformation']['pass'][:1] for path in metapaths]
        df['orbit'] = orbits
        # filter orbits
        if filter_orbit is not None:
            df = df[df.orbit == filter_orbit]

        df = df.set_index('date').sort_values('datetime')\
            [['datetime','orbit','mission','polarization','subswath','datapath','metapath','orbitpath','geometry']]

        err, warn = self.validate(df)
        #print ('err, warn', err, warn)
        assert not err, 'ERROR: Please fix all the issues listed above to continue'
        if warn:
            print ('NOTE: Please follow all the notes listed above')

        self.df = df
        # set first image as master
        self.master = self.df.index[0]
        
        self.set_dem(dem_filename)
        self.set_landmask(landmask_filename)
        
        # initialize empty pins list
        self.pins = []

    def validate(self, df=None):
        import numpy as np

        if df is None:
            df = self.df
        error = False
        warning = False

        # we can't merge together scenes from different missions
        missions = df.groupby('date')['mission'].unique().values
        missions = [len(mission) for mission in missions if len(mission)>1]
        if not len(missions) == 0:
            error = True
            print ('ERROR: Found multiple scenes for a single date from different missions')
        # note: df.unique() returns unsorted values so it would be 21 instead of expected 12
        subswaths = int(''.join(map(str,np.unique(df.subswath))))
        if not int(subswaths) in [1, 2, 3, 12, 23, 123]:
            error = True
            print (f'ERROR: Subswhats list {subswaths} incorrect. Allowed a single or sequential subswath numbers 1, 2, 3, 12, 23, 123')
        if not len(df.orbit.unique()) <= 1:
            error = True
            print ('ERROR: Only single orbit processing supported. Use any one ascending or descending')
        if not len(df.index.unique()) >= 2:
            error = True
            print ('ERROR: Two or more scenes required')
        daily_scenes = df.groupby(['date', 'subswath'])['datetime'].count().values.max()
        if daily_scenes > 1:
            warning = True
            print ('NOTE: Found multiple scenes for a single day, use function SBAS.reframe() to stitch the scenes')
        return error, warning   

    def backup(self, backup_dir, copy=False, debug=False):
        import os
        import shutil

        os.makedirs(backup_dir, exist_ok=True)

        # this optional file is dumped state, copy it if exists
        # auto-generated file can't be a symlink but user-defined symlink target should be copied
        filename = os.path.join(self.basedir, 'SBAS.pickle')
        if os.path.exists(filename):
            if debug:
                print ('DEBUG: copy', filename, backup_dir)
            shutil.copy2(filename, backup_dir, follow_symlinks=True)

        # these files required to continue the processing, do not remove and copy only
        filenames = [self.dem_filename, self.landmask_filename]
        for filename in filenames:
            # DEM and landmask can be not defined
            if filename is None:
                continue
            if debug:
                print ('DEBUG: copy', filename, backup_dir)
            shutil.copy2(filename, backup_dir, follow_symlinks=True)

        # these files are large and are not required to continue the processing
        filenames = []
        for record in self.df.itertuples():
            for filename in [record.datapath, record.metapath, record.orbitpath]:
                filenames += filename if isinstance(filename, list) else [filename]
        for filename in filenames:
            # orbit files can be not defined
            if filename is None:
                continue
            # copy and delete the original later to prevent cross-device links issues
            if debug:
                print ('DEBUG: copy', filename, backup_dir)
            shutil.copy2(filename, backup_dir, follow_symlinks=True)
            if not copy and self.basedir == os.path.dirname(filename):
                # when copy is not needed then delete files in work directory only
                if debug:
                    print ('DEBUG: remove', filename)
                os.remove(filename)

        if not copy:
            # mark as empty all the removed files
            for col in ['datapath','metapath','orbitpath']:
                self.df[col] = None
        return

    def set_master(self, master):
        """
        Set master image date
        """
        if not master in self.df.index:
            raise Exception('Master image not found')
        self.master = master
        return self

    def get_master(self, subswath=None):
        """
        Return dataframe master record(s) for all or only selected subswath
        """
        df = self.df.loc[[self.master]]
        if not subswath is None:
            df = df[df.subswath == subswath]
        assert len(df) > 0, f'Master record for subswath {subswath} not found'
        return df

    # enlist all the subswaths
    def get_subswaths(self):
        import numpy as np
        # note: df.unique() returns unsorted values so it would be 21 instead of expected 12
        return np.unique(self.df.subswath)
    
    def get_subswath(self, subswath=None):
        """
        Check and return subswath or return an unique subswath to functions which work with a single subswath only.
        """
        # detect all the subswaths
        subswaths = self.get_subswaths()
        assert subswath is None or subswath in subswaths, f'ERROR: subswath {subswath} not found'
        if subswath is not None:
            return subswath
        assert len(subswaths)==1, f'ERROR: multiple subswaths {subswaths} found, merge them first using SBAS.merge_parallel()'
        # define subswath
        return subswaths[0]
        
    # for precision orbit there is only single orbit per day
    # for approximate orbit 2 and maybe more orbits per day are possible
    # so check orbit file for for each subswath
    def download_orbits(self):
        from eof.download import download_eofs

        df = self.df[self.df['orbitpath'].isna()]
    
        # download all the misssed orbit files
        for record in df.itertuples():
            #print (date, mission)
            orbitpath = download_eofs([record.datetime], [record.mission], save_dir=self.basedir)[0]
            #print ('orbitpath', orbitpath)
            self.df.loc[self.df.datetime == record.datetime,'orbitpath'] = orbitpath

    def set_dem(self, dem_filename):
        import os
        if dem_filename is not None:
            self.dem_filename = os.path.relpath(dem_filename,'.')
        else:
            self.dem_filename = None
        return self

    def set_landmask(self, landmask_filename):
        import os
        if landmask_filename is not None:
            self.landmask_filename = os.path.relpath(landmask_filename,'.')
        else:
            self.landmask_filename = None
        return self

    # wrapper
    def download_dem(self, backend=None, **kwargs):
        if backend is None:
            return self.download_dem_old(**kwargs)
        elif backend == 'GMT':
            return self.download_dem_gmt(**kwargs)
        else:
            raise Exception(f'Unknown backend {backend}. Use None or GMT')
        

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

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small margin produces insufficient DEM not covers the defined area
    # https://docs.generic-mapping-tools.org/6.0/datasets/earth_relief.html
    def download_dem_gmt(self, product='SRTM1', resolution_meters=60, method='cubic', buffer_degrees=0.02, debug=False):
        """
        Use GMT server to download SRTM 1 or 3 arcsec data (@earth_relief_01s or @earth_relief_03s)
        Remove EGM96 geoid to make heights relative to WGS84
        Regrid to specified approximate resolution_meters (60m by default)
        """
        from pygmtsar import PRM
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

    def get_aligned(self, subswath=None, date=None):
        """
        Return dataframe aligned records (excluding master) for selected subswath
        """
        if date is None:
            idx = self.df.index.difference([self.master])
        else:
            idx = [date]
        df = self.df.loc[idx]
        if subswath is None:
            return df
        return df[df.subswath == subswath]

    def multistem_stem(self, subswath, dt=None):
        """
        Define stem and multistem using datetime    
        """
        from datetime import datetime

        # use master datetime if not defined
        if dt is None:
            dt = self.df.loc[self.master, 'datetime']

        stem = f'S1_{dt.strftime("%Y%m%d_%H%M%S")}_F{subswath}'
        multistem = f'S1_{dt.strftime("%Y%m%d")}_ALL_F{subswath}'
        return (multistem, stem)

    @staticmethod
    def annotation(filename):
        """
        Return XML scene annotation
        """
        import xmltodict

        with open(filename) as fd:
            # fix wrong XML tags to process cropped scenes
            # GMTSAR assemble_tops.c produces malformed xml
            # https://github.com/gmtsar/gmtsar/issues/354
            doc = xmltodict.parse(fd.read().replace('/></','></'))
        return doc

    def geoloc(self, filename=None):
        """
        Build approximate scene polygons using GCPs from XML scene annotation
        """
        from pygmtsar import PRM
        import numpy as np
        import pandas as pd
        import geopandas as gpd
        import os
    
        # for backward compatibility
        if filename is None:
            print ('NOTE: SBAS.geoloc() command is outdated. Use SBAS.to_dataframe().plot()')
            return pd.DataFrame({'longitude': [], 'latitude': [], 'pixel': ''})

        doc = self.annotation(filename)
        geoloc = doc['product']['geolocationGrid']['geolocationGridPointList']
        # check data consistency
        assert int(geoloc['@count']) == len(geoloc['geolocationGridPoint'])
        gcps = pd.DataFrame(geoloc['geolocationGridPoint']).applymap(lambda val : pd.to_numeric(val,errors='ignore'))
        # return approximate location as set of GCP
        return gpd.GeoDataFrame(gcps, geometry=gpd.points_from_xy(x=gcps.longitude, y=gcps.latitude))

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
        dem = xr.open_dataset(self.dem_filename, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
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

    def get_landmask(self, subswath=None, geoloc=False, buffer_degrees=0.02, inverse_geocode=False):
        import xarray as xr
        import os

        if self.landmask_filename is None:
            raise Exception('Set landmask first')

        # open DEM file and find the elevation variable
        # because sometimes grid includes 'crs' or other variables
        landmask = xr.open_dataset(self.landmask_filename, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
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
    
    def get_pins(self, subswath=None):
        """
        Return linear list of two pin coordinates for one or all subswaths. Use this list to easy plot the pins.
        """
    
        if subswath is None:
            # that's ok when pins are not defined here
            # all the pins useful for plotting only
            return self.pins
        else:
            # detect all the subswaths
            subswaths = self.get_subswaths()
            # check pins defined for all the subswaths
            assert len(self.pins)/4. == len(subswaths), f'ERROR: Pins are not defined for all the subswaths. \
                Found {len(self.pins)} pins for {subswaths} subswathss'
            assert subswath in subswaths, f'Subswath {subswath} not found'
            idx = subswaths.tolist().index(subswath)
            pins = self.pins[4*idx:4*(idx+1)]
            assert len(pins) == 4, f'ERROR: wrong number of pins detected. Found {len(pins)} for subswath {subswath}'
            return pins

    def set_pins(self, *args):
        """
        Estimate framed area using two pins on Sentinel-1 GCPs approximation.
        For each defined subswath the pins should be defined like to
            [x1, y1, x2, y2]
            [[x1, y1], [x2, y2]]
            [[x1, y1], None]
            [None, [x2, y2]]
            [None, None]
        The pins automatically ordering for ascending and descending orbits.
        """
        import numpy as np
        from shapely.geometry import Point
        # project pin location to boundary geometry
        def pip2pin(geom, lon, lat):
            pin = geom.boundary.interpolate(geom.boundary.project(Point(lon, lat)))
            return np.round(pin.x, 3), np.round(pin.y, 3)

        # detect all the subswaths
        subswaths = self.get_subswaths()
        if len(args) == 0:
            args = len(subswaths)*[None]
        # check input data to have a single argument for the each subswath
        assert len(args) == len(subswaths), f'Define pair of pins for the each subswath \
            {",".join(map(str,subswaths))} ({2*len(subswaths)} pins and {4*len(subswaths)} coordinates in total)'

        # iterate 1 to 3 subswaths
        allpins = []
        for (subswath, pins) in zip(subswaths, args):
            #print ('subswath, pins', subswath, pins)

            error = False
            warning = False

            if pins is None:
                pins = [None, None]
            if len(pins) == 4:
                pins = np.array(pins).reshape(2,2)
            assert len(pins) == 2, 'Define two pins as two pairs of lat,lon coordinates where pin2 is upper pin1 or use None'
            #print ('pins', pins)
            pin1 = pins[0]
            pin2 = pins[1]

            df = self.get_master(subswath)
            geom = df['geometry'].unary_union
            orbit = df['orbit'][0]

            # check the pins validity
            llmin, ltmin, llmax, ltmax = geom.envelope.bounds
            if not np.all(pin1) is None and not geom.intersects(Point(pin1[0], pin1[1])):
                print (f'ERROR subswath {subswath}: pin1 lays outside of master frame. Move the pin or set it to None and try again.')
                error = True
            if not np.all(pin2) is None and not geom.intersects(Point(pin2[0], pin2[1])):
                print (f'ERROR subswath {subswath}: pin2 lays outside of master frame. Move the pin or set it to None and try again.')
                error = True

            # check pin1
            if np.all(pin1) is None and orbit == 'A':
                # use right bottom corner
                print (f'NOTE subswath {subswath}: pin1 is not defined, master image corner coordinate will be used')
                warning = True
                #pin1 = [llmax, ltmin]
                pin1 = pip2pin(geom, llmax, ltmin)
            elif np.all(pin1) is None and orbit == 'D':
                # use right top corner
                print (f'NOTE subswath {subswath}: pin1 is not defined, master image corner coordinate will be used')
                warning = True
                #pin1 = [llmin, ltmin]
                pin1 = pip2pin(geom, llmin, ltmin)
            # check pin2
            if np.all(pin2) is None and orbit == 'A':
                # use left top corner
                print (f'NOTE subswath {subswath}: pin2 is not defined, master image corner coordinate will be used')
                warning = True
                #pin2 = [llmin, ltmax]
                pin2 = pip2pin(geom, llmin, ltmax)
            elif np.all(pin2) is None and orbit == 'D':
                # use left bottom corner
                print (f'NOTE subswath {subswath}: pin2 is not defined, master image corner coordinate will be used')
                warning = True
                #pin2 = [llmax, ltmax]
                pin2 = pip2pin(geom, llmax, ltmax)

            # check pins order
            if pin1[1] >= pin2[1]:
                print (f'ERROR subswath {subswath}: pin1 is upper than pin2. Fix to set pin1 at bottom and pin2 at top.')
                error = True
            # swap pins for Descending orbit
            if orbit == 'A':
                allpins += [pin1[0], pin1[1], pin2[0], pin2[1]]
            else:
                allpins += [pin2[0], pin2[1], pin1[0], pin1[1]]

        self.pins = allpins
        # pins are defined even is case of errors to have ability to plot them
        assert not error, 'ERROR: Please fix all the issues listed above to continue. Note: you are able to plot the pins.'

    def reframe(self, subswath, date, debug=False):
        """
        Estimate framed area using two pins using Sentinel-1 GCPs approximation.
        """
        from pygmtsar import PRM
        import numpy as np
        import shapely
        import os

        pins = self.get_pins(subswath)

        df = self.get_aligned(subswath, date)
        stem = self.multistem_stem(subswath, df['datetime'][0])[1]
        #print ('stem', stem)

        old_filename = os.path.join(self.basedir, f'{stem}')
        #print ('old_filename', old_filename)

        self.make_s1a_tops(subswath, date, debug=debug)

        prm = PRM.from_file(old_filename+'.PRM')
        tmpazi = prm.SAT_llt2rat([pins[0], pins[1], 0], precise=1)[1]
        if debug:
            print ('DEBUG: ','tmpazi', tmpazi)
        prm.shift_atime(tmpazi, inplace=True).update()
        azi1 = prm.SAT_llt2rat([pins[0], pins[1], 0], precise=1)[1] + tmpazi
        azi2 = prm.SAT_llt2rat([pins[2], pins[3], 0], precise=1)[1] + tmpazi
        if debug:
                print ('DEBUG: ','azi1', azi1, 'azi2', azi2)

        # Working on bursts covering $azi1 ($ll1) - $azi2 ($ll2)...
        #print ('assemble_tops', subswath, date, azi1, azi2, debug)
        self.assemble_tops(subswath, date, azi1, azi2, debug=debug)

        # Parse new .xml to define new scene name
        # like to 's1b-iw3-slc-vv-20171117t145922-20171117t145944-008323-00ebab-006'
        filename = os.path.splitext(os.path.split(df['datapath'][0])[-1])[0]
        head1 = filename[:15]
        tail1 = filename[-17:]
        xml_header = self.annotation(old_filename+'.xml')['product']['adsHeader']
        date_new = xml_header['startTime'][:10].replace('-','')
        t1 = xml_header['startTime'][11:19].replace(':','')
        t2 = xml_header['stopTime'][11:19].replace(':','')
        new_name = f'{head1}{date_new}t{t1}-{date_new}t{t2}-{tail1}'
        new_filename = os.path.join(self.basedir, new_name)
        #print ('new_filename', new_filename)

        # rename xml and tiff
        for ext in ['.tiff', '.xml']:
            if debug:
                print('DEBUG: rename', old_filename+ext, new_filename+ext)
            os.rename(old_filename+ext, new_filename+ext)

        # cleanup
        for fname in [old_filename+'.LED', old_filename+'.PRM']:
            if not os.path.exists(fname):
                continue
            if debug:
                print ('DEBUG: remove', fname)
            os.remove(fname)

        # update and return only one record
        df = df.head(1)
        df['datetime'] = self.text2date(f'{date_new}t{t1}', False)
        df['metapath'] = new_filename + '.xml'
        df['datapath'] = new_filename + '.tiff'
        # update approximate location
        gcps = self.geoloc(new_filename + '.xml').geometry
        df['geometry'] = shapely.geometry.MultiPoint(gcps).minimum_rotated_rectangle

        return df

    def reframe_parallel(self, dates=None, n_jobs=-1, **kwargs):
        from tqdm.auto import tqdm
        import joblib
        import pandas as pd

        if len(self.pins) == 0:
            # set auto-pins when the list is empty
            self.set_pins()

        if dates is None:
            dates = self.df.index.unique().values

        subswaths = self.get_subswaths()

        # process all the scenes
        with self.tqdm_joblib(tqdm(desc='Reframing', total=len(dates)*len(subswaths))) as progress_bar:
            records = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.reframe)(subswath, date, **kwargs) \
                                                     for date in dates for subswath in subswaths)

        self.df = pd.concat(records)
    
    def intf_ra2ll(self, subswath=None, grids=None, debug=False):
        """
        Geocoding function based on interferogram geocode matrix to call from open_grids(geocode=True)
        """
        from tqdm.auto import tqdm
        import joblib
        import xarray as xr
        import numpy as np
        import os

        # that's possible to miss the first argument subswath
        assert subswath is not None or grids is not None, 'ERROR: define input grids'
        if grids is None:
            grids = subswath
            subswath = None

        # helper check
        if not 'y' in grids.dims and 'x' in grids.dims:
            print ('NOTE: the grid is not in radar coordinates, miss geocoding')
            return grids

        # check if subswath exists or return a single subswath for None
        subswath = self.get_subswath(subswath)

        intf_ra2ll_file = os.path.join(self.basedir, f'F{subswath}_intf_ra2ll.grd')
        intf_ll2ra_file = os.path.join(self.basedir, f'F{subswath}_intf_ll2ra.grd')

        matrix_ra2ll = xr.open_dataarray(intf_ra2ll_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
        matrix_ll2ra = xr.open_dataarray(intf_ll2ra_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)

        # conversion works for a different 1st grid dimension size
        def ra2ll(grid):
            # input and transform grids should be the same
            grid = grid.reindex_like(matrix_ll2ra)
            # some interferograms have different y dimension and matrix has the largest
            # crop matrix y dimension when it is larger than current interferogram
            matrix_ra2ll_valid = xr.where(matrix_ra2ll<grid.size, matrix_ra2ll, -1)
            da_ll = xr.DataArray(np.where(matrix_ra2ll>=0, grid.values.reshape(-1)[matrix_ra2ll_valid], np.nan),
                coords=matrix_ra2ll_valid.coords)
            return da_ll

    #        def ra2ll(grid):
    #            return xr.DataArray(np.where(matrix_ra2ll>=0, grid.values.reshape(-1)[matrix_ra2ll], np.nan),
    #                coords=matrix_ra2ll.coords)

        # process single 2D raster
        if len(grids.dims) == 2:
            return ra2ll(grids)

        # process a set of 2D rasters
        with self.tqdm_joblib(tqdm(desc='Geocoding', total=len(grids))) as progress_bar:
            grids_ll = joblib.Parallel(n_jobs=-1)(joblib.delayed(ra2ll)(grids[item]) for item in range(len(grids)))
        grids_ll = xr.concat(grids_ll, dim=grids.dims[0])

        # add coordinates from original grids
        for coord in grids.coords:
            if coord in ['y', 'x']:
                continue
            grids_ll[coord] = grids[coord]

        return grids_ll

    def intf_ra2ll_matrix(self, subswath, intf_grids, debug=False):
        """
        Build interferogram geocoding matrix after interferogram processing required for open_grids(geocode=True)
        """
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np
        import os
    
        # use 2D grid grom the pairs stack
        # sometimes interferogram grids are different for one azimuth line so use the largest grid
        intf_grid = intf_grids.min('pair')

        trans_ra2ll_file = os.path.join(self.basedir, f'F{subswath}_trans_ra2ll.grd')
        intf_ra2ll_file  = os.path.join(self.basedir, f'F{subswath}_intf_ra2ll.grd')

        # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        trans = self.get_trans_dat(subswath)
        lon_min, lon_max = trans[:,3].min(),trans[:,3].max()
        lat_min, lat_max = trans[:,4].min(),trans[:,4].max()

        # read translation table for the full DEM area
        trans_ra2ll = xr.open_dataarray(trans_ra2ll_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)

        # build ra2ll translation matrix for interferogram coordinates and area only
        # each lan/lon cell has zero or one neighbour radar cell
        # each radar cell has one or multiple (on borders) neighbour lat/lon cells
        intf_ys, intf_xs = xr.broadcast(intf_grids.y, intf_grids.x)
        intf_yxs = np.stack([intf_ys.values.reshape(-1),intf_xs.values.reshape(-1)], axis=1)
        trans_yxs = np.stack([trans[:,1],trans[:,0]], axis=1)

        tree = cKDTree(intf_yxs, compact_nodes=False, balanced_tree=False)
        # use accurate distance limit as a half of the cell diagonal
        dy = intf_grids.y.diff('y')[0]
        dx = intf_grids.x.diff('x')[0]
        distance_limit = np.sqrt((dx/2.)**2 + (dy/2.)**2) + 1e-2
        d, inds = tree.query(trans_yxs, k = 1, distance_upper_bound=distance_limit, workers=8)

        # single integer index mask
        intf2trans = np.where(~np.isinf(d), inds, -1)
        # produce the same output array
        intf_ra2ll = xr.zeros_like(trans_ra2ll).rename('intf_ra2ll')
        intf_ra2ll.values = np.where(trans_ra2ll>=0, intf2trans[trans_ra2ll], -1)
        #assert intf_grid.size - 1 == intf_ra2ll.max(), 'ERROR: transform matrix and interferograms largest grid are different'
        assert intf_grid.size > intf_ra2ll.max(), \
            f'ERROR: transform matrix size {intf_grid.size} is too small for interferograms largest index {intf_ra2ll.max()}'
        # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
        intf_ra2ll.attrs['node_offset'] = 1
        # save to NetCDF file
        if os.path.exists(intf_ra2ll_file):
            os.remove(intf_ra2ll_file)
        intf_ra2ll.to_netcdf(intf_ra2ll_file, encoding={'intf_ra2ll': self.netcdf_compression}, engine=self.netcdf_engine)

    def ra2ll(self, subswath, debug=False):
        """
        Create radar to geographic coordinate transformation matrix for DEM grid using geocoding table trans.dat
        """
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np
        import os

        trans_ra2ll_file = os.path.join(self.basedir, f'F{subswath}_trans_ra2ll.grd')

        if os.path.exists(trans_ra2ll_file):
            os.remove(trans_ra2ll_file)

        # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        trans = self.get_trans_dat(subswath)
        lon_min, lon_max = trans[:,3].min(),trans[:,3].max()
        lat_min, lat_max = trans[:,4].min(),trans[:,4].max()

        #dem = xr.open_dataset(in_dem_gridfile)
        #dem = self.get_dem(geoloc=True)
        dem = self.get_dem(geoloc=True).sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))

        trans_latlons = np.stack([trans[:,4],trans[:,3]], axis=1)
        dem_lats, dem_lons = xr.broadcast(dem.lat,dem.lon)
        dem_latlons = np.stack([dem_lats.values.reshape(-1),dem_lons.values.reshape(-1)], axis=1)

        tree = cKDTree(trans_latlons, compact_nodes=False, balanced_tree=False)
        # use accurate distance limit as a half of the cell diagonal
        dlat = dem.lat.diff('lat')[0]
        dlon = dem.lon.diff('lon')[0]
        distance_limit = np.sqrt((dlat/2.)**2 + (dlon/2.)**2) + 1e-6
        d, inds = tree.query(dem_latlons, k = 1, distance_upper_bound=distance_limit, workers=8)

        # produce the same output array as dataset to be able to add global attributes
        trans_ra2ll = xr.zeros_like(dem).rename('trans_ra2ll')
        trans_ra2ll.values = np.where(~np.isinf(d), inds, -1).reshape(dem.shape)
        # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
        #trans_ra2ll.attrs['node_offset'] = 1
        # save to NetCDF file
        trans_ra2ll.to_netcdf(trans_ra2ll_file, encoding={'trans_ra2ll': self.netcdf_compression}, engine=self.netcdf_engine)

    # a single-step translation
    # see also two-step translation ra2ll & intf_ra2ll_matrix
    def intf_ll2ra_matrix(self, subswath, intf_grids, debug=False):
        """
        Create geographic to radar coordinate transformation matrix for DEM grid (F?_trans_ra2ll.grd)
        """
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np
        import os

        # use 2D grid grom the pairs stack
        # sometimes interferogram grids are different for one azimuth line so use the largest grid
        intf_grid = intf_grids.min('pair')

        trans_ra2ll_file = os.path.join(self.basedir, f'F{subswath}_trans_ra2ll.grd')
        intf_ll2ra_file  = os.path.join(self.basedir, f'F{subswath}_intf_ll2ra.grd')

        # to resolve opened NetCDF rewriting error
        if os.path.exists(intf_ll2ra_file):
            os.remove(intf_ll2ra_file)

        # read translation table for the full DEM area
        trans_ra2ll = xr.open_dataarray(trans_ra2ll_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)

        # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        trans = self.get_trans_dat(subswath)
        # convert trans lat, lon grid to the convertion table
        trans = trans[trans_ra2ll.values.reshape(-1),:]

        # the same indexing as in lat, lon conversion grid
        trans_yxs = np.stack([trans[:,1], trans[:,0]], axis=1)
        intf_ys, intf_xs = xr.broadcast(intf_grid.y, intf_grid.x)
        intf_yxs = np.stack([intf_ys.values.reshape(-1), intf_xs.values.reshape(-1)], axis=1)

        tree = cKDTree(trans_yxs, compact_nodes=False, balanced_tree=False)
        # use accurate distance limit as a half of the cell diagonal
        dy = intf_grid.y.diff('y')[0]
        dx = intf_grid.x.diff('x')[0]
        #distance_limit = np.sqrt((dy/2.)**2 + (dx/2.)**2) + 1e-2
        distance_limit = 100
        d, inds = tree.query(intf_yxs, k = 1, distance_upper_bound=distance_limit, workers=8)

        # produce the same output array as dataset to be able to add global attributes
        trans_ll2ra = xr.zeros_like(intf_grid).rename('intf_ll2ra')
        trans_ll2ra.values = np.where(~np.isinf(d), inds, -1).reshape(intf_grid.shape)

        if debug:
            # possible for merged subswaths due to subswaths offset
            undefined = (trans_ll2ra==-1).sum().item()
            print (f'DEBUG: inverse geocoding matrix has {undefined} undefined indices')
        # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
        trans_ll2ra.attrs['node_offset'] = 1
        # save to NetCDF file
        trans_ll2ra.to_netcdf(intf_ll2ra_file, encoding={'intf_ll2ra': self.netcdf_compression}, engine=self.netcdf_engine)

        return

    def intf_ll2ra(self, subswath=None, grids=None):
        """
        Inverse geocoding function based on interferogram geocode matrix to call from open_grids(inverse_geocode=True)
        """
        from tqdm.auto import tqdm
        import joblib
        import xarray as xr
        import numpy as np
        import os
        
        # that's possible to miss the first argument subswath
        assert subswath is not None or grids is not None, 'ERROR: define input grids'
        if grids is None:
            grids = subswath
            subswath = None
        
        # helper check
        if not 'lat' in grids.dims and 'lon' in grid.dims:
            print ('NOTE: the grid is not in geograpphic coordinates, miss geocoding')
            return grids

        # check if subswath exists or return a single subswath for None
        subswath = self.get_subswath(subswath)
            
        trans_ra2ll_file = os.path.join(self.basedir, f'F{subswath}_trans_ra2ll.grd')
        intf_ll2ra_file = os.path.join(self.basedir, f'F{subswath}_intf_ll2ra.grd')

        # read translation table for the full DEM area
        trans_ra2ll = xr.open_dataarray(trans_ra2ll_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
        # transform matrix
        matrix_ll2ra = xr.open_dataarray(intf_ll2ra_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)

        def ll2ra(grid):
            # transform input grid to the trans_ra2ll where the geocoding matrix defined
            # only nearest interpolation allowed to save values of binary masks
            return xr.DataArray(np.where(matrix_ll2ra>=0,
                                         grid.interp_like(trans_ra2ll, method='nearest').values.reshape(-1)[matrix_ll2ra],
                                         np.nan),
                coords=matrix_ll2ra.coords)

        # process single 2D raster
        if len(grids.dims) == 2:
            return ll2ra(grids)

        # process a set of 2D rasters
        with self.tqdm_joblib(tqdm(desc='Geocoding', total=len(grids))) as progress_bar:
            grids_ll = joblib.Parallel(n_jobs=-1)(joblib.delayed(ll2ra)(grids[item]) for item in range(len(grids)))
        grids_ll = xr.concat(grids_ll, dim=grids.dims[0])

        # add coordinates from original grids
        for coord in grids.coords:
            if coord in ['lat', 'lon']:
                continue
            grids_ll[coord] = grids[coord]

        return grids_ll

    def get_trans_dat(self, subswath=None):
        import numpy as np
        import os

        subswath = self.get_subswath(subswath)
        trans_dat_file  = os.path.join(self.basedir, f'F{subswath}_trans.dat')
        # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        trans_dat = np.fromfile(trans_dat_file, dtype=np.float64).reshape([-1,5])
        return trans_dat

    def trans_dat_parallel(self, n_jobs=-1, interactive=False, debug=False):
        import numpy as np
        import xarray as xr
        from tqdm.auto import tqdm
        import joblib
        import os

        # build trans.dat
        def SAT_llt2rat(subswath, ilat, ilon):
            # lazy dask array
            data = dem.data.blocks[ilat,ilon].compute()
            lats, lons = xr.broadcast(coordlat[ilat], coordlon[ilon])
            dem_data = np.column_stack([lons.values.ravel(), lats.values.ravel(), data.ravel()])
            return self.PRM(subswath).SAT_llt2rat(dem_data, precise=1, binary=True).reshape(-1,5)

        # process all the subswaths
        subswaths = self.get_subswaths()
        trans_dats = []
        for subswath in subswaths:
            trans_dat_file = os.path.join(self.basedir, f'F{subswath}_trans.dat')

            # cleanup before creating the new file
            if not interactive and os.path.exists(trans_dat_file):
                if debug:
                    print ('DEBUG: remove', trans_dat_file)
                os.remove(trans_dat_file)

            dem = self.get_dem(subswath, geoloc=True)
            lats, lons = dem.data.numblocks
            latchunks, lonchunks = dem.chunks
            coordlat = np.array_split(dem.lat, np.cumsum(latchunks))
            coordlon = np.array_split(dem.lon, np.cumsum(lonchunks))
            with self.tqdm_joblib(tqdm(desc=f'Radar Topography Computing', total=lats*lons)) as progress_bar:
                trans_dat_chunks = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(SAT_llt2rat)(subswath, ilat, ilon) \
                                               for ilat in range(lats) for ilon in range(lons))
            trans_dat = np.concatenate(trans_dat_chunks)
            if interactive:
                trans_dats.append(trans_dat)
            else:
                # save to binary file
                with open(trans_dat_file, 'wb') as f:
                    if debug:
                        print ('DEBUG: write', trans_dat_file)
                    f.write(trans_dat)
        if not interactive:
            return
        return trans_dats[0] if len(trans_dats)==1 else trans_dats

    # TODO: topo_ra.grd processing is not parallel but trans.dat computing only
    def topo_ra_parallel(self, idec=2, jdec=2, n_jobs=-1, debug=False):
        import numpy as np
        import xarray as xr
        from scipy.spatial import cKDTree
        from tqdm.auto import tqdm
        import joblib
        import os

        # build required trans.dat file
        self.trans_dat_parallel(n_jobs=n_jobs, debug=debug)

        # start a set of jobs together but not more than available cpu cores at once
        if n_jobs == -1:
            n_jobs = joblib.cpu_count()

        # as defined in dem2topo_ra.csh for Sentinel-1
        if debug:
            print (f'DEBUG: Topography range and azimuth decimation: {idec}/{jdec}')

        # process all the subswaths
        subswaths = self.get_subswaths()
        for subswath in subswaths:
            with self.tqdm_joblib(tqdm(desc=f'Radar Topography Gridding', total=1)) as progress_bar:
                topo_ra_file = os.path.join(self.basedir, f'F{subswath}_topo_ra.grd')
                # cleanup before creating the new file
                if os.path.exists(topo_ra_file):
                    if debug:
                        print ('DEBUG: remove', topo_ra_file)
                    os.remove(topo_ra_file)

                # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
                trans_dat = self.get_trans_dat(subswath)

                # build topo_ra
                XMAX, yvalid, num_patch = self.PRM(subswath).get('num_rng_bins', 'num_valid_az', 'num_patches')
                YMAX = yvalid * num_patch
                if debug:
                    print ('DEBUG: XMAX', XMAX, 'YMAX', YMAX)
                # use center pixel GMT registration mode
                rngs = np.arange(1, XMAX+1, idec)
                azis = np.arange(1, YMAX+1, jdec)
                grid_r, grid_a = np.meshgrid(rngs, azis)

                # find nearest DEM values for every radar grid point
                ras = np.stack([grid_r.reshape(-1), grid_a.reshape(-1)], axis=1)
                trans_ras = np.stack([trans_dat[:,0],trans_dat[:,1]], axis=1)
                tree = cKDTree(trans_ras, compact_nodes=False, balanced_tree=False)
                d, inds = tree.query(ras, k = 1, workers=n_jobs)
                if debug:
                    print (f'DEBUG: maximum distance is {np.max(d).round(1)} cells for subswath {subswath}')
                topo = trans_dat[:,2][inds].reshape(grid_r.shape)
                # flip vertically for GMTSAR compatibility reasons
                topo_ra = xr.DataArray(np.flipud(topo),
                                       dims=['y', 'x'],
                                       coords={'y': azis, 'x': rngs},
                                       name='z')
                # save to NetCDF file
                topo_ra.to_netcdf(topo_ra_file, encoding={'z': self.netcdf_compression}, engine=self.netcdf_engine)

                # update progress bar
                progress_bar.update(1)

    def get_topo_ra(self, subswath=None):
        import xarray as xr
        import dask.array
        import os

        if subswath is None:
            # process all the subswaths
            subswaths = self.get_subswaths()
        else:
            # process only one subswath
            subswaths = [subswath]

        topos = []
        for subswath in subswaths:
            topo_ra_file = os.path.join(self.basedir, f'F{subswath}_topo_ra.grd')

            # topography stored flipped vertically
            topo = xr.open_dataarray(topo_ra_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
            # flip vertically for GMTSAR compatibility reasons
            topo_fixed = xr.DataArray(dask.array.flipud(topo), coords=topo.coords, attrs=topo.attrs, name=topo.name)
            topos.append(topo_fixed)
            
        return topos[0] if len(topos)==1 else topos

    def sat_look_parallel(self, n_jobs=-1, interactive=False, debug=False):
        import numpy as np
        import xarray as xr
        from tqdm.auto import tqdm
        import joblib
        import os

        def SAT_look(subswath, ilat, ilon):
            # lazy dask arrays
            lats, lons = xr.broadcast(coordlat[ilat], coordlon[ilon])
            data = dem.sel(lat=xr.DataArray(lats.values.ravel()),
                                 lon=xr.DataArray(lons.values.ravel()),
                                 method='nearest').compute()
            coords = np.column_stack([lons.values.ravel(), lats.values.ravel(), data.values.ravel()])
            # look_E look_N look_U
            look = self.PRM(subswath).SAT_look(coords, binary=True)\
                                     .reshape(-1,6)[:,3:].astype(np.float32)
            # prepare output as xarray dataset
            dims = ['lat', 'lon']
            coords = coords={'lat': coordlat[ilat], 'lon':coordlon[ilon]}
            look_E = xr.DataArray(look[:,0].reshape(lats.shape), dims=dims, coords=coords, name='look_E')
            look_N = xr.DataArray(look[:,1].reshape(lats.shape), dims=dims, coords=coords, name='look_N')
            look_U = xr.DataArray(look[:,2].reshape(lats.shape), dims=dims, coords=coords, name='look_U')
            return xr.merge([look_E, look_N, look_U])

        # process all the subswaths
        subswaths = self.get_subswaths()
        dss = []
        for subswath in subswaths:
            intf_ra2ll_file = os.path.join(self.basedir, f'F{subswath}_intf_ra2ll.grd')
            sat_look_file   = os.path.join(self.basedir, f'F{subswath}_sat_look.grd')

            # cleanup before creating the new file
            if not interactive and os.path.exists(sat_look_file):
                if debug:
                    print ('DEBUG: remove', sat_look_file)
                os.remove(sat_look_file)

            dem = self.get_dem(subswath, geoloc=True)
            grid_ll = xr.open_dataarray(intf_ra2ll_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
            lats, lons = grid_ll.data.numblocks
            latchunks, lonchunks = grid_ll.chunks
            coordlat = np.array_split(grid_ll.lat, np.cumsum(latchunks))
            coordlon = np.array_split(grid_ll.lon, np.cumsum(lonchunks))
            with self.tqdm_joblib(tqdm(desc=f'SAT_look Computing', total=lats*lons)) as progress_bar:
                ds = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(SAT_look)(subswath, ilat, ilon) \
                                               for ilat in range(lats) for ilon in range(lons))
            # concatenate the chunks
            if xr.__version__ == '0.19.0':
                # for Google Colab
                ds = xr.merge(ds)
            else:
                # for modern xarray versions
                ds = xr.combine_by_coords(ds)

            if interactive:
                dss.append(ds)
            else:
                if debug:
                    print ('DEBUG: write NetCDF', sat_look_file)
                # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
                ds.attrs['node_offset'] = 1
                # save to NetCDF file
                ds.to_netcdf(sat_look_file,
                             encoding={var:self.netcdf_compression for var in ds.data_vars},
                             engine=self.netcdf_engine)
        if not interactive:
            return
        return dss[0] if len(dss)==1 else dss

    def get_sat_look(self, subswath=None):
        import xarray as xr
        import os

        if subswath is None:
            # process all the subswaths
            subswaths = self.get_subswaths()
        else:
            # process only one subswath
            subswaths = [subswath]

        sat_looks = []
        for subswath in subswaths:
            sat_look_file = os.path.join(self.basedir, f'F{subswath}_sat_look.grd')
            assert os.path.exists(sat_look_file), 'ERROR: satellite looks grid missed. Build it first using SBAS.sat_look_parallel()'
            sat_look = xr.open_dataset(sat_look_file, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
            sat_looks.append(sat_look)

        return sat_looks[0] if len(sat_looks)==1 else sat_looks

    def transforms(self, subswath=None, pairs=None, debug=False):
    
        assert pairs is not None or subswath is not None, 'ERROR: define pairs argument'
        if pairs is None and subswath is not None:
            pairs = subswath
            subswath = None
        
        subswath = self.get_subswath(subswath)
        if debug:
            print (f'DEBUG: build translation matrices for direct and inverse geocoding for subswath {subswath}')

        # build DEM grid coordinates transform matrix
        self.ra2ll(subswath, debug=debug)
    
        # transforms for interferogram grid
        grids = self.open_grids(pairs[:1], 'phasefilt')
        # build radar coordinates transformation matrix for the interferograms grid stack
        self.intf_ra2ll_matrix(subswath, grids, debug=debug)
        # build geographic coordinates transformation matrix for landmask and other grids
        self.intf_ll2ra_matrix(subswath, grids, debug=debug)
    
    # replacement for gmt grdfilter ../topo/dem.grd -D2 -Fg2 -I12s -Gflt.grd
    # use median decimation instead of average
    def get_topo_llt(self, subswath, degrees, geoloc=True, debug=False):
        import xarray as xr
        import numpy as np

        # add buffer around the cropped area for borders interpolation
        dem_area = self.get_dem(subswath, geoloc=geoloc)
        ny = int(np.round(degrees/dem_area.lat.diff('lat')[0]))
        nx = int(np.round(degrees/dem_area.lon.diff('lon')[0]))
        if debug:
            print ('DEBUG: DEM decimation','ny', ny, 'nx', nx)
        dem_area = dem_area.coarsen({'lat': ny, 'lon': nx}, boundary='pad').mean()

        lats, lons, z = xr.broadcast(dem_area.lat, dem_area.lon, dem_area)
        topo_llt = np.column_stack([lons.values.ravel(), lats.values.ravel(), z.values.ravel()])
        return topo_llt

    def offset2shift(self, xyz, rmax, amax, method='linear'):
        import xarray as xr
        import numpy as np
        from scipy.interpolate import griddata

        # use center pixel GMT registration mode
        rngs = np.arange(8/2, rmax+8/2, 8)
        azis = np.arange(4/2, amax+4/2, 4)
        grid_r, grid_a = np.meshgrid(rngs, azis)

        grid = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (grid_r, grid_a), method=method)
        da = xr.DataArray(np.flipud(grid), coords={'y': azis, 'x': rngs}, name='z')
        return da

    def to_dataframe(self):
        return self.df

    def assemble_tops(self, subswath, date, azi_1, azi_2, debug=False):
        """
        Usage: assemble_tops azi_1 azi_2 name_stem1 name_stem2 ... output_stem

        Example: assemble_tops 1685 9732 s1a-iw1-slc-vv-20150706t135900-20150706t135925-006691-008f28-001
            s1a-iw1-slc-vv-20150706t135925-20150706t135950-006691-008f28-001
            s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001

        Output:s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001.xml
            s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001.tiff

        Note: output files are bursts that covers area between azi_1 and azi_2, set them to 0s to output all bursts
        """
        import numpy as np
        import os
        import subprocess

        df = self.get_aligned(subswath, date)
        #print ('scenes', len(df))

        # assemble_tops requires the same path to xml and tiff files
        datadirs = [os.path.split(path)[:-1] for path in df['datapath']]
        metadirs = [os.path.split(path)[:-1] for path in df['metapath']]
        if not datadirs == metadirs:
            # in case when the files placed in different directories we need to create symlinks for them
            datapaths = []
            for datapath, metapath in zip(df['datapath'], df['metapath']):
                for filepath in [datapath, metapath]:
                    filename = os.path.split(filepath)[-1]
                    relname = os.path.join(self.basedir, filename)
                    if os.path.exists(relname) or os.path.islink(relname):
                        os.remove(relname)
                    os.symlink(os.path.relpath(filepath, self.basedir), relname)
                datapaths.append(os.path.splitext(filename)[0])
        else:
            datapaths = [os.path.relpath(path, self.basedir)[:-5] for path in df['datapath']]
        #print ('datapaths', datapaths)
        stem = self.multistem_stem(subswath, df['datetime'][0])[1]

        # round values and convert to strings
        azi_1 = np.round(azi_1).astype(int).astype(str)
        azi_2 = np.round(azi_2).astype(int).astype(str)

        argv = ['assemble_tops', azi_1, azi_2] + datapaths + [stem]
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: assemble_tops', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: assemble_tops', stdout_data.decode('ascii'))

        return

    def ext_orb_s1a(self, subswath, stem, date=None, debug=False):
        import os
        import subprocess

        if date is None or date == self.master:
            df = self.get_master(subswath)
        else:
            df = self.get_aligned(subswath, date)

        orbit = os.path.relpath(df['orbitpath'][0], self.basedir)

        argv = ['ext_orb_s1a', f'{stem}.PRM', orbit, stem]
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: ext_orb_s1a', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: ext_orb_s1a', stdout_data.decode('ascii'))

        return
    
    # produce LED and PRM in basedir
    # when date=None work on master image
    def make_s1a_tops(self, subswath, date=None, mode=0, rshift_fromfile=None, ashift_fromfile=None, debug=False):
        """
        Usage: make_slc_s1a_tops xml_file tiff_file output mode dr.grd da.grd
         xml_file    - name of xml file 
         tiff_file   - name of tiff file 
         output      - stem name of output files .PRM, .LED, .SLC 
         mode        - (0) no SLC; (1) center SLC; (2) high SLCH and lowSLCL; (3) output ramp phase
        """
        import os
        import subprocess

        #or date == self.master
        if date is None:
            df = self.get_master(subswath)
            # for master image mode should be 1
            mode = 1
        else:
            df = self.get_aligned(subswath, date)

        # TODO: use subswath
        xmlfile = os.path.relpath(df['metapath'][0], self.basedir)
        datafile = os.path.relpath(df['datapath'][0], self.basedir)
        stem = self.multistem_stem(subswath, df['datetime'][0])[1]

        argv = ['make_s1a_tops', xmlfile, datafile, stem, str(mode)]
        if rshift_fromfile is not None:
            argv.append(rshift_fromfile)
        if ashift_fromfile is not None:
            argv.append(ashift_fromfile)
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: make_s1a_tops', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: make_s1a_tops', stdout_data.decode('ascii'))

        self.ext_orb_s1a(subswath, stem, date, debug=debug)

        return

    # aligning for master image
    def stack_ref(self, subswath, debug=False):
        import xarray as xr
        import numpy as np
        import os
        from pygmtsar import PRM

    #        err, warn = self.validate()
    #        #print ('err, warn', err, warn)
    #        assert not err and not warn, 'ERROR: Please fix all the issues listed above to continue'

        master_line = list(self.get_master(subswath).itertuples())[0]
        #print (master_line)

        # for master image
        multistem, stem = self.multistem_stem(subswath, master_line.datetime)
        path_stem = os.path.join(self.basedir, stem)
        path_multistem = os.path.join(self.basedir, multistem)

        # generate PRM, LED, SLC
        self.make_s1a_tops(subswath, debug=debug)

        PRM.from_file(path_stem + '.PRM')\
            .set(input_file = path_multistem + '.raw')\
            .update(path_multistem + '.PRM', safe=True)

        self.ext_orb_s1a(subswath, multistem, debug=debug)

        # recalculate after ext_orb_s1a
        earth_radius = PRM.from_file(path_multistem + '.PRM')\
            .calc_dop_orb(inplace=True).update().get('earth_radius')

    # aligning for secondary image
    def stack_rep(self, subswath, date=None, degrees=12.0/3600, debug=False):
        import xarray as xr
        import numpy as np
        import os
        from pygmtsar import PRM
        
        # temporary filenames to be removed
        cleanup = []

        master_line = list(self.get_master(subswath).itertuples())[0]
        multistem, stem = self.multistem_stem(subswath, master_line.datetime)
        #print (master_line)

        # define master image parameters
        master = self.PRM(subswath).sel('earth_radius').set(stem=stem, multistem=multistem)

        # prepare coarse DEM for alignment
        # 12 arc seconds resolution is enough, for SRTM 90m decimation is 4x4
        topo_llt = self.get_topo_llt(subswath, degrees=degrees)
        #topo_llt.shape

        line = list(self.get_aligned(subswath, date).itertuples())[0]
        multistem, stem = self.multistem_stem(subswath, line.datetime)
        #print (line)

        # define relative filenames for PRM
        stem_prm    = os.path.join(self.basedir, stem + '.PRM')
        mstem_prm   = os.path.join(self.basedir, multistem + '.PRM')
        master_prm  = os.path.join(self.basedir, master.get("stem") + '.PRM')
        mmaster_prm = os.path.join(self.basedir, master.get("multistem") + '.PRM')

        # TODO: define 1st image for line, in the example we have no more
        tmp_da = 0

        # generate PRM, LED
        self.make_s1a_tops(subswath, date, debug=debug)

        # compute the time difference between first frame and the rest frames
        t1, prf = PRM.from_file(stem_prm).get('clock_start', 'PRF')
        t2      = PRM.from_file(stem_prm).get('clock_start')
        nl = int((t2 - t1)*prf*86400.0+0.2)
        #echo "Shifting the master PRM by $nl lines..."

        # Shifting the master PRM by $nl lines...
        # shift the super-masters PRM based on $nl so SAT_llt2rat gives precise estimate
        prm1 = PRM.from_file(master_prm)
        prm1.set(prm1.sel('clock_start' ,'clock_stop', 'SC_clock_start', 'SC_clock_stop') + nl/prf/86400.0)
        tmp_prm = prm1

        # compute whether there are any image offset
        #if tmp_da == 0:
        # tmp_prm defined above from {master}.PRM
        prm1 = tmp_prm.calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug)
        prm2 = PRM.from_file(stem_prm).calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug).update()
        lontie,lattie = prm1.SAT_baseline(prm2, debug=debug).get('lon_tie_point', 'lat_tie_point')
        tmp_am = prm1.SAT_llt2rat(coords=[lontie, lattie, 0], precise=1, debug=debug)[1]
        tmp_as = prm2.SAT_llt2rat(coords=[lontie, lattie, 0], precise=1, debug=debug)[1]
        # bursts look equal to rounded result int(np.round(...))
        tmp_da = int(tmp_as - tmp_am)
        #print ('tmp_am', tmp_am, 'tmp_as', tmp_as, 'tmp_da', tmp_da)

        # in case the images are offset by more than a burst, shift the super-master's PRM again
        # so SAT_llt2rat gives precise estimate
        if abs(tmp_da) >= 1000:
            prf = tmp_prm.get('PRF')
            tmp_prm.set(tmp_prm.sel('clock_start' ,'clock_stop', 'SC_clock_start', 'SC_clock_stop') - tmp_da/prf/86400.0)
            #raise Exception('TODO: Modifying master PRM by $tmp_da lines...')

        # tmp.PRM defined above from {master}.PRM
        prm1 = tmp_prm.calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug)
        tmpm_dat = prm1.SAT_llt2rat(coords=topo_llt, precise=1, debug=debug)
        prm2 = PRM.from_file(stem_prm).calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug)
        tmp1_dat = prm2.SAT_llt2rat(coords=topo_llt, precise=1, debug=debug)

        # get r, dr, a, da, SNR table to be used by fitoffset.csh
        offset_dat0 = np.hstack([tmpm_dat, tmp1_dat])
        func = lambda row: [row[0],row[5]-row[0],row[1],row[6]-row[1],100]
        offset_dat = np.apply_along_axis(func, 1, offset_dat0)

        # define radar coordinates extent
        rmax, amax = PRM.from_file(stem_prm).get('num_rng_bins','num_lines')

        # prepare the offset parameters for the stitched image
        # set the exact borders in radar coordinates
        par_tmp = offset_dat[(offset_dat[:,0]>0) & (offset_dat[:,0]<rmax) & (offset_dat[:,2]>0) & (offset_dat[:,2]<amax)]
        par_tmp[:,2] += nl
        if abs(tmp_da) >= 1000:
            par_tmp[:,2] -= tmp_da
            par_tmp[:,3] += tmp_da

        # prepare the rshift and ashift look up table to be used by make_s1a_tops
        # use tmp_dat instead of offset_dat
        r_xyz = offset_dat[:,[0,2,1]]
        a_xyz = offset_dat[:,[0,2,3]]

        r_grd = self.offset2shift(r_xyz, rmax, amax)
        r_grd_filename = stem_prm[:-4]+'_r.grd'
        r_grd.to_netcdf(r_grd_filename)
        # drop the temporary file at the end of the function
        cleanup.append(r_grd_filename)

        a_grd = self.offset2shift(a_xyz, rmax, amax)
        a_grd_filename = stem_prm[:-4]+'_a.grd'
        a_grd.to_netcdf(a_grd_filename)
        # drop the temporary file at the end of the function
        cleanup.append(a_grd_filename)

        # generate the image with point-by-point shifts
        # note: it removes calc_dop_orb parameters from PRM file
        # generate PRM, LED
        self.make_s1a_tops(subswath,
                           date=line.Index, mode=1,
                           rshift_fromfile=f'{stem}_r.grd',
                           ashift_fromfile=f'{stem}_a.grd',
                           debug=debug)

        # need to update shift parameter so stitch_tops will know how to stitch
        PRM.from_file(stem_prm).set(PRM.fitoffset(3, 3, offset_dat)).update()

        # echo stitch images together and get the precise orbit
        # use stitch_tops tmp.stitchlist $stem to merge images

        # the raw file does not exist but it works
        PRM.from_file(stem_prm)\
            .set(input_file = f'{multistem}.raw')\
            .update(mstem_prm, safe=True)

        self.ext_orb_s1a(subswath, multistem, date=line.Index, debug=debug)

        # Restoring $tmp_da lines shift to the image... 
        PRM.from_file(mstem_prm).set(ashift=0 if abs(tmp_da) < 1000 else tmp_da, rshift=0).update()

        # that is safe to rewrite source files
        prm1 = PRM.from_file(mmaster_prm)
        prm1.resamp(PRM.from_file(mstem_prm),
                    alignedSLC_tofile=mstem_prm[:-4]+'.SLC',
                    interp=1, debug=debug
        ).to_file(mstem_prm)

        PRM.from_file(mstem_prm).set(PRM.fitoffset(3, 3, par_tmp)).update()

        PRM.from_file(mstem_prm).calc_dop_orb(master.get('earth_radius'), 0, inplace=True, debug=debug).update()
        
        # cleanup
        for filename in cleanup:
            #if os.path.exists(filename):
            os.remove(filename)

    def stack_parallel(self, dates=None, n_jobs=-1, **kwargs):
        from tqdm.auto import tqdm
        import joblib

        if dates is None:
            dates = list(self.get_aligned().index.unique())

        subswaths = self.get_subswaths()

        # prepare master image
        #self.stack_ref()
        with self.tqdm_joblib(tqdm(desc='Reference', total=len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.stack_ref)(subswath, **kwargs) for subswath in subswaths)

        # prepare secondary images
        with self.tqdm_joblib(tqdm(desc='Aligning', total=len(dates)*len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.stack_rep)(subswath, date, **kwargs) \
                                           for date in dates for subswath in subswaths)

    def intf(self, subswath, pair, **kwargs):
        from pygmtsar import PRM
        import os

        # extract dates from pair
        date1, date2 = pair

        prm_ref = self.PRM(subswath, date1)
        prm_rep = self.PRM(subswath, date2)

        topo_ra_file = os.path.join(self.basedir, f'F{subswath}_topo_ra.grd')
        #print ('SBAS intf kwargs', kwargs)
        prm_ref.intf(prm_rep,
                     basedir=self.basedir,
                     topo_ra_fromfile = topo_ra_file,
                     **kwargs)

    def intf_parallel(self, pairs, n_jobs=-1, **kwargs):
        import pandas as pd
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        from joblib.externals import loky
        import os

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        subswaths = self.get_subswaths()

        # this way does not work properly for long interferogram series
        #with self.tqdm_joblib(tqdm(desc='Interferograms', total=len(pairs))) as progress_bar:
        #    joblib.Parallel(n_jobs=-1)(joblib.delayed(self.intf)(pair, **kwargs) for pair in pairs)

        # start a set of jobs together but not more than available cpu cores at once
        if n_jobs == -1:
            n_jobs = joblib.cpu_count()
        n_chunks = int(np.ceil(len(pairs)/n_jobs))
        chunks = np.array_split(pairs, n_chunks)
        #print ('n_jobs', n_jobs, 'n_chunks', n_chunks, 'chunks', [len(chunk) for chunk in chunks])
        with tqdm(desc='Interferograms', total=len(pairs)*len(subswaths)) as pbar:
            for chunk in chunks:
                loky.get_reusable_executor(kill_workers=True).shutdown(wait=True)
                with joblib.parallel_backend('loky', n_jobs=n_jobs, inner_max_num_threads=1):
                    joblib.Parallel()(joblib.delayed(self.intf)(subswath, pair, **kwargs) \
                        for subswath in subswaths for pair in chunk)
                pbar.update(len(chunk)*len(subswaths))

        # backward compatibility wrapper
        # for a single subswath don't need to call SBAS.merge_parallel()
        # for subswaths merging and total coordinate transformation matrices creation 
        if len(subswaths) == 1:
            # build geo transform matrices for interferograms
            self.transforms(subswaths[0], pairs)
        
    # stem_tofile + '.PRM' generating
    def merge_swath(self, conf, grid_tofile, stem_tofile, debug=False):
        import subprocess

        argv = ['merge_swath', '/dev/stdin', grid_tofile, stem_tofile]
        if debug:
            print ('DEBUG: argv', argv)
            print ('DEBUG: conf:', f'\n{conf}')
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='ascii')
        stdout_data, stderr_data = p.communicate(input=conf)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: merge_swath', stderr_data)
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: merge_swath', stdout_data)

        return

    def merge(self, pair, grid, debug=False):
        from pygmtsar import PRM
        import os

        record2multistem = lambda record: self.multistem_stem(record.subswath, record.datetime)[0]
        fullname = lambda filename: os.path.join(self.basedir, filename)

        # extract dates from pair
        date1, date2 = pair
        #print (date1, date2)

        # records should be sorted by datetime that's equal to sorting by date and subswath
        multistems1 = self.get_aligned(None, date1).apply(record2multistem, axis=1)
        multistems2 = self.get_aligned(None, date2).apply(record2multistem, axis=1)
        if len(multistems1) == 1:
            # only one subswath found, merging is not possible
            return

        config = []
        subswaths = []
        cleanup = []
        for (multistem1, multistem2) in zip(multistems1, multistems2):
            # F2_20220702_20220714_phasefilt.grd
            prm1_filename = fullname(multistem1 + '.PRM')
            prm2_filename = fullname(multistem1 + '.PRM')
            prm_filename  = fullname(multistem1 + f'_{grid}.PRM')

            prm1 = PRM.from_file(prm1_filename)
            prm2 = PRM.from_file(prm2_filename)
            rshift = prm2.get('rshift')
            #print ('rshift', rshift)
            #assert rshift == 0, 'rshift is not equal to zero for master PRM'
            fs1 = prm1.get('first_sample')
            fs2 = prm2.get('first_sample')
            #print ('fs1, fs2', fs1, fs2)
            #assert fs1 == fs2, 'first_sample is not equal for master and repeat PRM'
            prm = prm1.set(rshift=rshift, first_sample=fs2 if fs2 > fs1 else fs1).to_file(prm_filename)

            subswath = int(multistem1[-1:])
            subswaths.append(subswath)
            dt1 = multistem1[3:11]
            dt2 = multistem2[3:11]
            #print (multistem1, multistem2, fullname(multistem1))
            grid_fromfile = fullname(f'F{subswath}_{dt1}_{dt2}_{grid}.grd')
            cleanup.append(grid_fromfile)
            #print (prm_filename, grid_filename)
            config.append(':'.join([prm_filename, grid_fromfile]))
        config = '\n'.join(config)
        subswaths = int(''.join(map(str,subswaths)))
        # F23_20220702_20220714_phasefilt.grd
        grid_tofile = fullname(f'F{subswaths}_{dt1}_{dt2}_{grid}.grd')
        # F23_20220702_20220714_phasefilt
        tmp_stem_tofile = fullname(f'F{subswaths}_{dt1}_{dt2}_{grid}')
        #print ('grid_tofile', grid_tofile)
        #print (config)
        # S1_20220702_ALL_F23 without extension
        stem_tofile = fullname(f'S1_{dt1}_ALL_F{subswaths}')

        # use temporary well-qualified stem file name to prevent parallel processing conflicts
        self.merge_swath(config, grid_tofile, tmp_stem_tofile, debug=debug)
        # different pairs and grids generate the same PRM file, replace it silently
        if debug:
            print ('DEBUG: replace', f'{tmp_stem_tofile}.PRM', f'{stem_tofile}.PRM')
        os.replace(f'{tmp_stem_tofile}.PRM', f'{stem_tofile}.PRM')

        # cleanup - files should exists as these are processed above
        for filename in cleanup:
            if debug:
                print ('DEBUG: remove', filename)
            os.remove(filename)

    def merge_parallel(self, pairs, grids = ['phasefilt', 'corr'], n_jobs=-1, **kwargs):
        from tqdm.auto import tqdm
        import joblib
        import pandas as pd
        import geopandas as gpd

        # merging is not applicable to a single subswath
        # for this case coordinate transformation matrices already built in SBAS.intf_parallel()
        subswaths = self.get_subswaths()
        if len(subswaths) == 1:
            return
        
        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        with self.tqdm_joblib(tqdm(desc=f'Merging Subswaths', total=len(pairs)*len(grids))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.merge)(pair, grid, **kwargs) \
                                           for pair in pairs for grid in grids)

        df = self.df.groupby(self.df.index).agg({'datetime': min, 'orbit': min, 'mission': min, 'polarization':min,
                                            'subswath': lambda s: int(''.join(map(str,list(s)))),
                                            'datapath': lambda p: list(p),
                                            'metapath': lambda p: list(p),
                                            'orbitpath': min,
                                            'geometry': lambda g: g.unary_union
                                           })
        self.df = gpd.GeoDataFrame(df)
    
        # build topo_ra and geo transform matrices for the merged interferograms
        self.topo_ra_parallel()
        self.transforms(pairs)

    def baseline_table(self):
        import pandas as pd
        import numpy as np

        # use any subswath (actually, the 1st one) to produce the table
        subswath = self.get_subswaths()[0]
        # select unique dates to process multiple subswaths
        dates = np.unique(self.df.index)

        # after merging use unmerged subswath PRM files
        prm_ref = self.PRM(subswath, singleswath=True)
        data = []
        for date in dates:
            # after merging use unmerged subswath PRM files
            prm_rep = self.PRM(subswath, date, singleswath=True)
            ST0 = prm_rep.get('SC_clock_start')
            DAY = int(ST0 % 1000)
            YR = int(ST0/1000) - 2014
            YDAY = YR * 365 + DAY
            #print (f'YR={YR}, DAY={DAY}')
            BPL, BPR = prm_ref.SAT_baseline(prm_rep).get('B_parallel', 'B_perpendicular')
            data.append({'date':date, 'ST0':ST0, 'YDAY':YDAY, 'BPL':BPL, 'BPR':BPR})
            #print (date, ST0, YDAY, BPL, BPR)
        return pd.DataFrame(data).set_index('date')

    # returns sorted baseline pairs
    def baseline_pairs(self, days=100, meters=150, invert=False):
        import numpy as np
        import pandas as pd
     
        tbl = self.baseline_table()
        data = []
        for line1 in tbl.itertuples():
        #for line1 in tbl.loc[['2015-01-21']].itertuples():
            for line2 in tbl.itertuples():
            #for line2 in tbl.loc[['2015-03-10']].itertuples():
                #print (line1, line2)
                if not (line1.YDAY < line2.YDAY and line2.YDAY - line1.YDAY < days):
                    continue
                if not (abs(line1.BPR - line2.BPR)< meters):
                    continue

                if not invert:
                    data.append({'ref_date':line1.Index, 'rep_date': line2.Index,
                                 'ref_timeline': np.round(line1.YDAY/365.25+2014, 2), 'ref_baseline': np.round(line1.BPR, 2),
                                 'rep_timeline': np.round(line2.YDAY/365.25+2014, 2), 'rep_baseline': np.round(line2.BPR, 2)})
                else:
                    data.append({'ref_date':line2.Index, 'rep_date': line1.Index,
                                 'ref_timeline': np.round(line2.YDAY/365.25+2014, 2), 'ref_baseline': np.round(line2.BPR, 2),
                                 'rep_timeline': np.round(line1.YDAY/365.25+2014, 2), 'rep_baseline': np.round(line1.BPR, 2)})

        return pd.DataFrame(data).sort_values(['ref_date', 'rep_date'])

    @staticmethod
    def triplets2pairs(triplets, pairs):
        import pandas as pd
    
        data = []
        for triplet in triplets.itertuples():
            data.append({'ref_date': triplet.A, 'rep_date': triplet.B})
            data.append({'ref_date': triplet.B, 'rep_date': triplet.C})
            data.append({'ref_date': triplet.A, 'rep_date': triplet.C})
        tripairs = pd.DataFrame(data).sort_values(['ref_date', 'rep_date']).drop_duplicates()
        idx = tripairs.set_index(['ref_date', 'rep_date']).index
        return pairs.set_index(['ref_date', 'rep_date']).loc[idx].reset_index()

    # returns sorted baseline triplets
    @staticmethod
    def pairs2triplets(pairs, invert=False):
        import pandas as pd

        data = []
        pairs_a = pairs
        for line_a in pairs_a.itertuples():
            #print (line_a)
            date_a_ref = line_a.ref_date
            date_a_rep = line_a.rep_date
            pairs_b = pairs[pairs.ref_date==date_a_rep]
            for line_b in pairs_b.itertuples():
                #print (line_b)
                date_b_ref = line_b.ref_date
                date_b_rep = line_b.rep_date
                pairs_c = pairs[(pairs.rep_date==date_b_rep)&(pairs.ref_date==date_a_ref)]
                for line_c in pairs_c.itertuples():
                    #print (line_c)
                    date_c_ref = line_c.ref_date
                    date_c_rep = line_c.rep_date
                    #print (date_a_ref, date_a_rep, date_b_rep)
                    data.append({'A': date_a_ref, 'B': date_a_rep, 'C': date_b_rep})
        return pd.DataFrame(data).sort_values(['A', 'B', 'C'])

    def PRM(self, subswath=None, date=None, multi=True, singleswath=False):
        """
        multi=True/False - open multistem or stem file
        singleswath=False/True - open a single-digit subswath instead of a multi-digit (merged) one
            single-digit subswath exists always while multi-digit exists only for interferogram pair references
        """
        from pygmtsar import PRM
        import os

        # check if subswath exists or return a single subswath for None
        subswath = self.get_subswath(subswath)

        if date is None or date == self.master:
            line = self.get_master(subswath)
        else:
            line = self.get_aligned(subswath, date)
        #print (line)
        # to build sbas table and pairs after merging use unmerged subswath PRM files
        if singleswath and len(str(subswath))>1:
            subswath = int(str(subswath)[0])
        multistem, stem = self.multistem_stem(subswath, line.datetime[0])
        if multi:
            stem = multistem
        filename = os.path.join(self.basedir, f'{stem}.PRM')
        #print (filename)
        return PRM.from_file(filename)

    def unwrap_parallel(self, pairs, n_jobs=-1, **kwargs):
        import pandas as pd
        from tqdm.auto import tqdm
        import joblib
        import os

        def unwrap(subswath, pair, **kwargs):
            # define unique tiledir name for parallel processing
            if 'conf' in kwargs:
                dirname = f'F{subswath}_{"_".join(pair).replace("-","")}_snaphu_tiledir'
                dirpath = os.path.join(self.basedir, dirname)
                kwargs['conf'] += f'    TILEDIR {dirpath}'
            return self.unwrap(subswath, pair, **kwargs)

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        subswaths = self.get_subswaths()

        with self.tqdm_joblib(tqdm(desc='Unwrapping', total=len(pairs)*len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(unwrap)(subswath, pair, interactive=False, **kwargs) \
                                           for subswath in subswaths for pair in pairs)

    # -s for SMOOTH mode and -d for DEFO mode when DEFOMAX_CYCLE should be defined in the configuration
    # DEFO mode (-d) and DEFOMAX_CYCLE=0 is equal to SMOOTH mode (-s)
    # https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu_man1.html
    def unwrap(self, subswath, pair, threshold=None, conf=None, func=None, mask=None, conncomp=False,
               interactive=True, debug=False):
        import xarray as xr
        import numpy as np
        import os
        import subprocess

        if threshold is None:
            # set to very low value but still exclude 0 (masked areas)
            threshold = 1e-6

        # convert user-defined mask to binary mask (NaN values converted to 0)
        if mask is not None:
            assert self.is_ra(mask), 'ERROR: mask should be defined in radar coordinates'
            # crop mask to minimum extent
            binmask = self.cropna(xr.where(mask>0, 1, np.nan)).fillna(0).astype(bool)
        else:
            binmask = 1

        if conf is None:
            conf = self.PRM(subswath).snaphu_config()

        # extract dates from pair
        date1, date2 = pair

        basename = os.path.join(self.basedir, f'F{subswath}_{date1}_{date2}_').replace('-','')
        #print ('basename', basename)

        # input data grids
        phase_filename = basename + 'phasefilt.grd'
        corr_filename = basename + 'corr.grd'
        # output data grids
        unwrap_filename = basename + 'unwrap.grd'
        conncomp_filename = basename + 'conncomp.grd'

        # SNAPHU input files
        phase_in = basename + 'unwrap.phase'
        corr_in = basename + 'unwrap.corr'
        # SNAPHU output files
        unwrap_out = basename + 'unwrap.out'
        conncomp_out = basename + 'conncomp.out'

        phase = xr.open_dataarray(phase_filename, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
        corr = xr.open_dataarray(corr_filename, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
        if mask is not None:
            phase = phase.reindex_like(binmask)
            corr = corr.reindex_like(binmask)

        # cleanup from previous runs
        for tmp_file in [phase_in, corr_in, unwrap_out, conncomp_out,
                         conncomp_filename, unwrap_filename]:
            #print ('tmp_file', tmp_file)
            if os.path.exists(tmp_file):
                if debug:
                    print ('DEBUG: remove', tmp_file)
                os.remove(tmp_file)

        # prepare SNAPHU input files
        # NaN values are not allowed for SNAPHU phase input file
        phase.where(~np.isnan(phase),0).astype(np.float32).values.tofile(phase_in)
        # apply threshold and binary mask
        corrmasked = corr.where(binmask & (corr>=threshold))
        # NaN values are not allowed for SNAPHU correlation input file
        corrmasked.fillna(0).astype(np.float32).values.tofile(corr_in)

        # launch SNAPHU binary (NaNs are not allowed for input but returned in output)
        argv = ['snaphu', phase_in, str(phase.shape[1]), '-c', corr_in,
                '-f', '/dev/stdin', '-o', unwrap_out, '-d']
        if conncomp:
            argv.append('-g')
            argv.append(conncomp_out)
        if debug:
            argv.append('-v')
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='ascii', bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=conf)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: snaphu', stderr_data)
        if debug:
            print ('DEBUG: snaphu', stdout_data)

        if conncomp:
            # convert to grid the connected components from SNAPHU output as is (UCHAR)
            values = np.fromfile(conncomp_out, dtype=np.ubyte).reshape(phase.shape)
            conn = xr.DataArray(values, phase.coords, name='conncomp')

        # convert to grid unwrapped phase from SNAPHU output applying postprocessing
        values = np.fromfile(unwrap_out, dtype=np.float32).reshape(phase.shape)
        #values = np.frombuffer(stdout_data, dtype=np.float32).reshape(phase.shape)
        unwrap = xr.DataArray(values, phase.coords)
        # apply user-defined function for post-processing
        if func is not None:
            unwrap = func(corrmasked, unwrap)
        # apply binary mask after the post-processing to completely exclude masked regions
        # NaN values allowed in the output grid, assign userfriendly name for the output grid
        unwrap = unwrap.where(binmask>0).rename('phase')

        for tmp_file in [phase_in, corr_in, unwrap_out, conncomp_out]:
            if os.path.exists(tmp_file):
                if debug:
                    print ('DEBUG: remove', tmp_file)
                os.remove(tmp_file)

        if interactive:
            return (unwrap, conn) if conncomp else unwrap

        # not interactive mode, save all the results to disk
        if conncomp:
            conn.to_netcdf(conncomp_filename, encoding={'conncomp': self.netcdf_compression}, engine=self.netcdf_engine)
        # save to NetCDF file
        unwrap.to_netcdf(unwrap_filename, encoding={'phase': self.netcdf_compression}, engine=self.netcdf_engine)

    #gmt grdmath unwrap_mask.grd $wavel MUL -79.58 MUL = los.grd
    def los_displacement_mm(self, unwraps):
        # constant is negative to make LOS = -1 * range change
        # constant is (1000 mm) / (4 * pi)
        scale = -79.58 * self.PRM().get('radar_wavelength')
        los_disp = scale*unwraps
        return los_disp

    def incidence_angle(self, subswath=None):
        import xarray as xr
        import numpy as np

        subswath = self.get_subswath(subswath)
        sat_look = self.get_sat_look(subswath)
        incidence_ll = np.arctan2(np.sqrt(sat_look.look_E**2 + sat_look.look_N**2), sat_look.look_U).rename('incidence_angle')
        return incidence_ll

    def vertical_displacement_mm(self, unwraps):
        import numpy as np
    
        assert self.is_geo(unwraps), 'ERROR: unwrapped phase defined in radar coordinates'
        
        los_disp = self.los_displacement_mm(unwraps)
        incidence_ll = self.incidence_angle()
        return los_disp/np.cos(incidence_ll)

    def eastwest_displacement_mm(self, unwraps):
        import numpy as np
    
        # this displacement is not symmetrical for the orbits due to scene geometries
        orbit = self.df.orbit.unique()[0]
        sign = 1 if orbit == 'D' else -1
        los_disp = self.los_displacement_mm(unwraps)
        incidence_ll = self.incidence_angle()
        return sign * los_disp/np.sin(incidence_ll)

    # TODO: use the function for open_grids() and do the same for open_grid() too
    def get_filenames(self, subswath, pairs, name, add_subswath=True):
        import pandas as pd
        import numpy as np
        import os

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values
        else:
            pairs = np.asarray(pairs)

        if add_subswath == True:
            prefix = f'F{subswath}_'
        else:
            prefix = ''

        filenames = []
        if len(pairs.shape) == 1:
            for date in sorted(pairs):
                filename = os.path.join(self.basedir, f'{prefix}{name}_{date}.grd'.replace('-',''))
                filenames.append(filename)
        elif len(pairs.shape) == 2:
            for pair in pairs:
                filename = os.path.join(self.basedir, f'{prefix}{pair[0]}_{pair[1]}_{name}.grd'.replace('-',''))
                filenames.append(filename)
        return filenames
    
    # returns all grids in basedir by mask or grids by dates and name
    # Backward-compatible open_grids() returns list of grids fot the name or a single grid for a single subswath
    def open_grids(self, pairs, name, geocode=False, inverse_geocode=False,  mask=None, func=None,
                   crop_valid=False, add_subswath=True, chunks=None, n_jobs=-1):
        import pandas as pd
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        import os

        assert not(geocode and inverse_geocode), 'ERROR: Only single geocoding option can be applied'

        if chunks is None:
            chunks = self.netcdf_chunksize

        # iterate all the subswaths
        subswaths = self.get_subswaths()

        # for backward compatibility
        if isinstance(mask, bool):
            print ('NOTE: mask argument changed from boolean to dataarray for SBAS.open_grid() function call')
            mask = None
        if mask is not None:
            assert len(subswaths) == 1, 'ERROR: mask can be applied to a single subswath only'
            nanmask = xr.where((mask == 0)|(np.isnan(mask)), np.nan, 1)

        if pairs is None:
            # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
            pass
        elif isinstance(pairs, pd.DataFrame):
            pairs = pairs.values
        else:
            pairs = np.asarray(pairs)

        def postprocess(da, subswath):
            if self.is_ra(da) and geocode:
                da = self.intf_ra2ll(subswath, da)
            elif self.is_geo(da) and inverse_geocode:
                da = self.intf_ll2ra(subswath, da)
            if func is not None:
                if isinstance(func, list):
                    for f in func:
                        da = f(da)
                else:
                    da = func(da)
            if mask is not None:
                assert self.is_same(mask, da), 'ERROR: mask defined in different coordinates'
                da = nanmask * da
            return da
    
        dass = []
        for subswath in subswaths:
            if add_subswath == True:
                prefix = f'F{subswath}_'
            else:
                prefix = ''

            das = []
            if pairs is None:
                # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
                filename = os.path.join(self.basedir, f'{prefix}{name}.grd')
                #print ('filename', filename)
                da = xr.open_dataarray(filename, engine=self.netcdf_engine, chunks=chunks)
                das  = postprocess(da, subswath)
            elif len(pairs.shape) == 1:
                # read all the grids from files
                for date in sorted(pairs):
                    filename = os.path.join(self.basedir, f'{prefix}{name}_{date}.grd'.replace('-',''))
                    #print (date, filename)
                    da = xr.open_dataarray(filename, engine=self.netcdf_engine, chunks=chunks)
                    das.append(da)

                # post-processing on a set of 2D rasters
                with self.tqdm_joblib(tqdm(desc='Loading', total=len(das))) as progress_bar:
                    das = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(postprocess)(da, subswath) for da in das)

                # prepare to stacking
                das = [da.expand_dims('date') for da in das]

                # 2nd dimension sizes must be the same
                #sizes = len(np.unique([da[0].shape[1] for da in das]))
                #assert sizes==1, f'Dates grids have different {da.dims[0]} dimensions'
                # 1st dimension can be different a bit, use minimal size
                #sizes = np.unique([da[0].shape[0] for da in das])
                #das = xr.concat(das, dim='date')[:,:sizes[0],:]

                # allow stack to be extended to largest 1st dimension size
                # to be sure all code work well for this case
                # so user is able to load grids by his own way
                das = xr.concat(das, dim='date')
                das['date'] = sorted(pairs)
            elif len(pairs.shape) == 2:
                # read all the grids from files
                for pair in pairs:
                    filename = os.path.join(self.basedir, f'{prefix}{pair[0]}_{pair[1]}_{name}.grd'.replace('-',''))
                    #print (filename)
                    da = xr.open_dataarray(filename, engine=self.netcdf_engine, chunks=chunks)
                    das.append(da)

                # post-processing on a set of 2D rasters
                with self.tqdm_joblib(tqdm(desc='Loading', total=len(das))) as progress_bar:
                    das = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(postprocess)(da, subswath) for da in das)

                # prepare to stacking
                das = [da.expand_dims('pair') for da in das]

                # 2nd dimension sizes must be the same
                #sizes = len(np.unique([da[0].shape[1] for da in das]))
                #assert sizes==1, f'Pairs grids have different {da.dims[0]} dimensions'

                # 1st dimension can be different a bit, use minimal size
                # this works but it is not robust if user loads grids by other way
                #sizes = np.unique([da[0].shape[0] for da in das])
                #das = xr.concat(das, dim='pair')[:,:sizes[0],:]

                # allow stack to be extended to largest 1st dimension size
                # to be sure all code work well for this case
                # so user is able to load grids by his own way
                das = xr.concat(das, dim='pair')
                das['pair'] = [f'{pair[0]} {pair[1]}' for pair in pairs]
                das['ref']  = xr.DataArray([pair[0] for pair in pairs], dims='pair')
                das['rep']  = xr.DataArray([pair[1] for pair in pairs], dims='pair')
            else:
                raise Exception('Use single or two columns Pandas dataset or array as "pairs" argument')

            # crop NaNs
            if crop_valid:
                das = self.cropna(das)
            dass.append(das)

        return dass[0] if len(dass) == 1 else dass

    def pixel_size(self, grid=(1, 4), average=True):
        import xarray as xr
        import numpy as np

        outs = []
        for subswath in self.get_subswaths():
            # pixel size in meters
            azi_px_size, rng_px_size = self.PRM(subswath).pixel_size()
            # raster pixels decimation
            if isinstance(grid, xr.DataArray):
                dy = grid.y.diff('y')[0].item()
                dx = grid.x.diff('x')[0].item()
            else:
                dy, dx = grid
            outs.append((np.round(azi_px_size*dy,1), np.round(rng_px_size*dx,1)))
        if average:
            pxs = np.asarray(outs)
            return (pxs[:,0].mean(), pxs[:,1].mean())
        else:
            return outs[0] if len(outs) == 1 else outs

    #decimator = lambda da: da.coarsen({'y': 2, 'x': 2}, boundary='trim').mean()
    def pixel_decimator(self, resolution_meters=60, debug=False):
        import numpy as np
    
        dy, dx = self.pixel_size()
        yy, xx = int(np.round(resolution_meters/dy)), int(np.round(resolution_meters/dx))
        if debug:
            print (f'DEBUG: average per subswaths ground pixel size in meters: y={dy}, x={dx}')
        if debug:
            print (f"DEBUG: decimator = lambda da: da.coarsen({{'y': {yy}, 'x': {xx}}}, boundary='trim').mean()")
        return lambda da: da.coarsen({'y': yy, 'x': xx}, boundary='trim').mean()

    def detrend_parallel(self, pairs, n_jobs=-1, **kwargs):
        from tqdm.auto import tqdm
        import joblib
        import pandas as pd

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        subswaths = self.get_subswaths()

        # process all the scenes
        with self.tqdm_joblib(tqdm(desc='Detrending', total=len(pairs)*len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.detrend)(subswath, pair, interactive=False, **kwargs) \
                                                     for subswath in subswaths for pair in pairs)

    def detrend(self, subswath, pair=None, wavelength=None, topo_ra=None, truncate=3.0, approximate=True,
                fit_intercept=True, fit_dem=True, fit_coords=True, interactive=True, debug=False):
        """
        Detrend and gaussian filtering on unwrapped interferogram in radar coordinates, see for details
            https://github.com/gmtsar/gmtsar/issues/98
            https://github.com/gmtsar/gmtsar/issues/411
        TODO: decimate grid for long-wavelength filtering for faster calculation.
        """
        import xarray as xr
        import numpy as np
        from sklearn.linear_model import LinearRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        import os

        if pair is None:
            pair = subswath
            subswath = self.get_subswath()
    
        # extract dates from pair
        date1, date2 = pair
        basename = os.path.join(self.basedir, f'F{subswath}_{date1}_{date2}_').replace('-','')

        # input data grid
        phase_filename = basename + 'unwrap.grd'
        # output data grid
        detrend_filename = basename + 'detrend.grd'
    
        #print ('basename', basename)
        phase = xr.open_dataarray(phase_filename, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)

        # topography grid defined in radar coordinates for larger area and different spacing
        if topo_ra is None:
            topo_ra = self.get_topo_ra()

        def postprocessing(out):
            #print ('interactive', interactive)
            if interactive:
                return out.astype(np.float32).rename('phase')
            # save to NetCDF file
            if os.path.exists(detrend_filename):
                os.remove(detrend_filename)
            out.astype(np.float32).rename('phase').to_netcdf(detrend_filename, encoding={'phase': self.netcdf_compression}, engine=self.netcdf_engine)

        # raster pixel spacing
        dy, dx = self.pixel_size(phase)
        ty, tx = self.pixel_size(topo_ra)

        # unify topo grid
        dec_topo_y = int(np.round(dy/ty))
        dec_topo_x = int(np.round(dx/tx))
        if debug:
            print ('DEBUG: Decimate topography to data grid dec_y, dec_x', dec_topo_y, dec_topo_x)
        # use xr.zeros_like to prevent the target grid coordinates modifying
        topo = topo_ra\
            .coarsen({'y': dec_topo_y, 'x': dec_topo_x}, boundary='pad').mean()\
            .reindex_like(xr.zeros_like(phase), method='nearest')

        if wavelength is not None:
            if debug:
                print ('DEBUG: Gaussian filtering for the both data and topography')
            phase.values -= self.gaussian(phase, wavelength, truncate=truncate, approximate=approximate, debug=debug)
            topo.values  -= self.gaussian(topo,  wavelength, truncate=truncate, approximate=approximate, debug=debug)

        # prepare raster
        y = phase.values.reshape(-1)
        nanmask = np.isnan(y)
        # prepare regression variable
        Y = y[~nanmask]
        if fit_coords or fit_dem:
            # prepare coordinates for X regression variable
            yy, xx = xr.broadcast(phase.y, phase.x)
            ys = yy.values.reshape(-1)
            xs = xx.values.reshape(-1)

        if fit_dem:
            # prepare topography for X regression variable
            zs = topo.values.reshape(-1)
            zys = zs*ys
            zxs = zs*xs

        if fit_dem and fit_coords:
            if debug:
                print ('DEBUG: Detrend topography and coordinates')
            X = np.column_stack([zys[~nanmask], zxs[~nanmask], ys[~nanmask], xs[~nanmask], zs[~nanmask]])
        elif fit_dem:
            if debug:
                print ('DEBUG: Detrend topography only')
            X = np.column_stack([zys[~nanmask], zxs[~nanmask], zs[~nanmask]])
        elif fit_coords:
            if debug:
                print ('DEBUG: Detrend coordinates only')
            X = np.column_stack([ys[~nanmask], xs[~nanmask]])
        elif fit_intercept:
            if debug:
                print ('DEBUG: Remove mean value only')
            return postprocessing(phase - phase.mean())
        else:
            if debug:
                print ('DEBUG: No detrending')
            return postprocessing(phase)

        # build prediction model with or without plane removal (fit_intercept)
        regr = make_pipeline(StandardScaler(), LinearRegression(fit_intercept=fit_intercept))
        regr.fit(X, Y)
        #model = np.nan * xr.zeros_like(da)
        # even for input dask array return numpy array
        model = xr.DataArray(np.nan * np.zeros(phase.shape), coords=phase.coords)
        model.values.reshape(-1)[~nanmask] = regr.predict(X)

        return postprocessing(phase - model)

    def make_gaussian_filter(self, range_dec, azi_dec, wavelength, debug=False):
        """
        Wrapper for PRM.make_gaussian_filter() and sonamed command line tool. Added for development purposes only.
        """
        import numpy as np

        gauss_dec, gauss_string = self.PRM().make_gaussian_filter(range_dec, azi_dec, wavelength, debug=debug)
        coeffs = [item for line in gauss_string.split('\n') for item in line.split('\t') if item != '']
        # x,y dims order
        shape = np.array(coeffs[0].split(' ')).astype(int)
        # y,x dims order
        matrix = np.array(coeffs[1:]).astype(float).reshape((shape[1],shape[0]))
        return (gauss_dec, matrix)

    def sbas_parallel(self, pairs, mask=None, detrended=True, data_stack=None, corr_stack=None, n_jobs=-1):
        import xarray as xr
        import numpy as np
        import pandas as pd
        from tqdm.auto import tqdm
        import joblib
        import os

        # compress 3d output following the processing blocks
        netcdf_compression = self.netcdf_compression.copy()
        netcdf_compression['chunksizes'] = (1, self.netcdf_chunksize, self.netcdf_chunksize)

        model_filename = os.path.join(self.basedir, 'disp.grd')
    
        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values
        else:
            pairs = np.asarray(pairs)
        # define all the dates as unique reference and repeat dates
        dates = np.unique(pairs.flatten())
    
        # source grids lazy loading
        if corr_stack is None:
            corr_stack = self.open_grids(pairs, 'corr')
        if data_stack is None:
            gridname = 'detrend' if detrended else 'unwrap'
            data_stack = self.open_grids(pairs, gridname)

        # crop correlation grid like to unwrap grid which may be defined smaller
        corr_stack = corr_stack.reindex_like(data_stack)
        
        # mask can be sparse and limit work area
        if mask is not None:
            data_stack = xr.where(mask>0, data_stack.reindex_like(mask), np.nan)
            corr_stack   = xr.where(mask>0, corr_stack.reindex_like(mask),   np.nan)
    
        # here are one row for every interferogram and one column for every date
        matrix = []
        for pair in pairs:
            mrow = [date>pair[0] and date<=pair[1] for date in dates]
            matrix.append(mrow)
        matrix = np.stack(matrix).astype(int)

        # single-pixel processing function
        def fit(x, w, matrix):
            #return np.zeros(5)
            # ignore pixels where correlation is not defined
            if np.any(np.isnan(w)):
                return np.nan * np.zeros(matrix.shape[1])
            # fill nans as zeroes and set corresponding weight to 0
            nanmask = np.where(np.isnan(x))
            if nanmask[0].size > 0:
                # function arguments are read-only
                x = x.copy()
                w = w.copy()
                x[nanmask] = 0
                w[nanmask] = 0
                # check if x has enough valid values
                if x.size - nanmask[0].size < matrix.shape[1]:
                    return np.nan * np.zeros(matrix.shape[1])
            # least squares solution
            W = (w/np.sqrt(1-w**2))
            model = np.linalg.lstsq(matrix * W[:,np.newaxis], x * W, rcond=None)
            #print ('model', model)
            return model[0]

        # xarray wrapper
        models = xr.apply_ufunc(
            fit,
            data_stack.chunk(dict(pair=-1)),
            corr_stack.chunk(dict(pair=-1)),
            input_core_dims=[['pair'],['pair']],
            exclude_dims=set(('pair',)),
            dask='parallelized',
            vectorize=True,
            output_dtypes=[np.float32],
            output_core_dims=[['date']],
            dask_gufunc_kwargs={'output_sizes': {'date': len(dates)}},
            kwargs={'matrix': matrix}
        )
        # define dates axis
        models['date'] = dates
        # set the stack index to be first
        models = models.transpose('date',...)
        # cleanup
        if os.path.exists(model_filename):
            os.remove(model_filename)
    
        ts, ys, xs = models.data.blocks.shape
        assert ts == 1, 'Date chunks count should be equal to 1'
        tchunks, ychunks, xchunks = models.chunks
        coordt = models.date
        coordy = np.array_split(models.y, np.cumsum(ychunks))
        coordx = np.array_split(models.x, np.cumsum(xchunks))
    
        def func(iy, ix):
            chunk_filename = os.path.join(self.basedir, f'disp_chunk_{iy}_{ix}.grd')
            if os.path.exists(chunk_filename):
                os.remove(chunk_filename)
            # lazy dask array
            data = models.data.blocks[0,iy,ix]
            # wrap the array
            da = xr.DataArray(data,
                              dims=['date','y','x'],
                              coords={'date': coordt, 'y':coordy[iy], 'x':coordx[ix]})\
                 .rename('displacement')
            # compute and save to NetCDF using chunks of original coordinates
            da.to_netcdf(chunk_filename,
                         unlimited_dims=['y','x'],
                         encoding={'displacement': netcdf_compression},
                         engine=self.netcdf_engine,
                         compute=True)
            return chunk_filename
    
        # process all the chunks
        with self.tqdm_joblib(tqdm(desc='Computing', total=ys*xs)) as progress_bar:
            filenames = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(func)(iy, ix) \
                                                     for iy in range(ys) for ix in range(xs))
    
        # rebuild the datasets to user-friendly format
        das = [xr.open_dataarray(f, engine=self.netcdf_engine, chunks=self.netcdf_chunksize) for f in filenames]
        if xr.__version__ == '0.19.0':
            # for Google Colab
            das = xr.merge(das)
        else:
            # for modern xarray versions
            das = xr.combine_by_coords(das)

        # add subswath prefix
        subswath = self.get_subswath()

        def output(dt):
            filename = os.path.join(self.basedir, f'F{subswath}_disp_{dt}.grd'.replace('-',''))
            if os.path.exists(filename):
                os.remove(filename)
            das.sel(date=dt).to_netcdf(filename,
                        encoding={'displacement': self.netcdf_compression},
                        engine=self.netcdf_engine)

        # saving all the grids
        with self.tqdm_joblib(tqdm(desc='Saving', total=len(das.date))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(output)(dt) for dt in das.date.values)

        # cleanup
        for filename in filenames:
            if os.path.exists(filename):
                os.remove(filename)

    #intf.tab format:   unwrap.grd  corr.grd  ref_id  rep_id  B_perp 
    def intftab(self, baseline_pairs):
        import numpy as np
        import datetime

        # return a single subswath for None
        subswath = self.get_subswath()

        outs = []
        for line in baseline_pairs.itertuples():
            #print (line)
            ref = line.ref_date.replace('-','')
            jref = datetime.datetime.strptime(line.ref_date, '%Y-%m-%d').strftime('%Y%j')
            rep = line.rep_date.replace('-','')
            jrep = datetime.datetime.strptime(line.rep_date, '%Y-%m-%d').strftime('%Y%j')
            bperp = np.round(line.rep_baseline - line.ref_baseline, 2)
            outs.append(f'F{subswath}_{ref}_{rep}_unwrap.grd F{subswath}_{ref}_{rep}_corr.grd {jref} {jrep} {bperp}')
        return '\n'.join(outs) + '\n'

    def scenetab(self, baseline_pairs):
        import numpy as np
        import datetime

        mst = datetime.datetime.strptime(self.master, '%Y-%m-%d').strftime('%Y%j')
        #print (self.master, mst)
        outs = []
        for line in baseline_pairs.itertuples():
            #print (line)
            ref = datetime.datetime.strptime(line.ref_date, '%Y-%m-%d').strftime('%Y%j')
            yday_ref = np.round((line.ref_timeline - 2014)*365.25)
            outs.append(f'{ref} {yday_ref}')
            rep = datetime.datetime.strptime(line.rep_date, '%Y-%m-%d').strftime('%Y%j')
            yday_rep = np.round((line.rep_timeline - 2014)*365.25)
            outs.append(f'{rep} {yday_rep}')
        outs = np.unique(outs)
        return '\n'.join([out for out in outs if out.split(' ')[0]==mst]) + '\n' + \
               '\n'.join([out for out in outs if out.split(' ')[0]!=mst]) + '\n'

    def sbas(self, baseline_pairs, smooth=0, atm=0, debug=False):
        """
         USAGE: sbas intf.tab scene.tab N S xdim ydim [-atm ni] [-smooth sf] [-wavelength wl] [-incidence inc] [-range -rng] [-rms] [-dem]

         input:
          intf.tab             --  list of unwrapped (filtered) interferograms:
           format:   unwrap.grd  corr.grd  ref_id  rep_id  B_perp
          scene.tab            --  list of the SAR scenes in chronological order
           format:   scene_id   number_of_days
           note:     the number_of_days is relative to a reference date
          N                    --  number of the interferograms
          S                    --  number of the SAR scenes
          xdim and ydim        --  dimension of the interferograms
          -smooth sf           --  smoothing factors, default=0
          -atm ni              --  number of iterations for atmospheric correction, default=0(skip atm correction)
          -wavelength wl       --  wavelength of the radar wave (m) default=0.236
          -incidence theta     --  incidence angle of the radar wave (degree) default=37
          -range rng           --  range distance from the radar to the center of the interferogram (m) default=866000
          -rms                 --  output RMS of the data misfit grids (mm): rms.grd
          -dem                 --  output DEM error (m): dem.grd

         output:
          disp_##.grd          --  cumulative displacement time series (mm) grids
          vel.grd              --  mean velocity (mm/yr) grids

         example:
          sbas intf.tab scene.tab 88 28 700 1000
        """
        import os
        import subprocess
        import numpy as np
        import math
        import glob
        import datetime

        # cleanup
        for filename in glob.glob(os.path.join(self.basedir, 'disp*.grd')):
            os.remove(filename)
        filename = os.path.join(self.basedir, 'vel.grd')
        if os.path.exists(filename):
            os.remove(filename)

        unwrap = self.open_grids(baseline_pairs[['ref_date', 'rep_date']][:1], 'unwrap')[0]
        dem = self.get_dem(geoloc=True)
        prm = self.PRM()

        #N=$(wc -l intf.in   | cut -d ' ' -f1)
        #S=$(wc -l scene.tab | cut -d ' ' -f1)

        N = len(baseline_pairs)
        S = len(np.unique(list(baseline_pairs['ref_date']) + list(baseline_pairs['rep_date'])))

        #bounds = self.geoloc().dissolve().envelope.bounds.values[0]
        llmin, ltmin, llmax, ltmax = self.get_master().dissolve().envelope.bounds.values[0].round(3)
        lon0 = (llmin + llmax)/2
        lat0 = (ltmin + ltmax)/2
        elevation0 = float(dem.sel(lat=lat0, lon=lon0, method='nearest'))
        #print ('coords',lon0, lat0, elevation0)
        _,_,_,look_E,look_N,look_U = prm.SAT_look([lon0, lat0, elevation0])
        #print ('satlook', _,_,_,look_E,look_N,look_U)
        incidence = math.atan2(math.sqrt(float(look_E)**2 + float(look_N)**2), float(look_U))*180/np.pi

        ydim, xdim = unwrap.shape

        xmin = int(unwrap.x.min())
        xmax = int(unwrap.x.max())
        near_range, rng_samp_rate, wavelength = prm.get('near_range', 'rng_samp_rate', 'radar_wavelength')
        # calculation below requires bc utility
        rng_pixel_size = 300000000 / rng_samp_rate / 2
        rng = np.round(rng_pixel_size * (xmin+xmax) /2 + near_range)

        intf_tab = self.intftab(baseline_pairs)
        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(intf_tab, 'ascii'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        scene_tab = self.scenetab(baseline_pairs)
        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(scene_tab, 'ascii'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        argv = ['sbas', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}',
                str(N), str(S), str(xdim), str(ydim), '-atm', str(atm), '-smooth', str(smooth),
                '-wavelength', str(wavelength), '-incidence', str(incidence), '-range', str(rng),
                '-rms', '-dem']
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=self.basedir, encoding='ascii')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: sbas', stderr_data)
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: sbas', stdout_data)

        # fix output grid filenames
        for date in np.unique(np.concatenate([baseline_pairs.ref_date,baseline_pairs.rep_date])):
            jdate = datetime.datetime.strptime(date, '%Y-%m-%d').strftime('%Y%j')
            date = date.replace('-','')
            filename1 = os.path.join(self.basedir, f'disp_{jdate}.grd')
            filename2 = os.path.join(self.basedir, f'disp_{date}.grd')
            if os.path.exists(filename1):
                if debug:
                    print ('DEBUG: rename', filename1, filename2)
                os.rename(filename1, filename2)
            #print (jdate, date)

        return
