#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to scan *xml files and orbits and make data.in file like to prep_data_linux.csh & prep_data.csh tools
#import pytest

class SBAS:

    # save to NetCDF grid
    compression = dict(zlib=True, complevel=3, chunksizes=[128,128])

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

    @staticmethod
    def nearest_grid(in_grid, search_radius_pixels=300):
        from PRM import PRM
        return PRM.nearest_grid(in_grid, search_radius_pixels)

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

    def __repr__(self):
        return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    def __init__(self, datadir, dem_filename=None, basedir='.',
                filter_orbit=None, filter_mission=None, filter_subswath=None, filter_polarization=None):
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
        assert filter_subswath is None or filter_subswath in [1,2,3], \
            'ERROR: use number 1 or 2 or 3 for subswath filter'
        assert filter_polarization is None or filter_polarization in ['VV','VH','HH','HV'], \
            'ERROR: use VV or VH or HH or HV for polarization filter'

        if dem_filename is None:
            self.dem_filename = None
        else:
            self.dem_filename = os.path.relpath(dem_filename, '.')
        #print ('dem_filename', self.dem_filename)

        # processing directory
        if basedir is None:
            self.basedir = '.'
        else:
            # (re)create basedir
            if os.path.exists(basedir):
                shutil.rmtree(basedir)
            os.makedirs(basedir, exist_ok=True)
            self.basedir = basedir

        if filter_polarization is None:
            filter_polarization = '??'
        if filter_subswath is None:
            filter_subswath  = '?'
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

        # filter mission
        if filter_mission is not None:
            orbit_path_pattern = f'{filter_mission.upper()}_OPER_AUX_*.EOF'
        else:
            orbit_path_pattern = 'S1?_OPER_AUX_*.EOF'
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
        df['subswath'] = [filename.split('-')[1] for filename in datanames]
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

    def validate(self, df=None):
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
        if not len(df.subswath.unique()) <= 1:
            error = True
            print ('ERROR: Only single subswath processing supported. Use any one iw1, iw2, or iw3')
        if not len(df.orbit.unique()) <= 1:
            error = True
            print ('ERROR: Only single orbit processing supported. Use any one ascending or descending')
        if not len(df.index.unique()) >= 2:
            error = True
            print ('ERROR: Two or more scenes required')
        daily_scenes = df.groupby('date')['datetime'].count().values.max()
        if daily_scenes > 1:
            warning = True
            print ('NOTE: Found multiple scenes for a single day, use function SBAS.reframe() to stitch the scenes')
        return error, warning
    
    def download_orbits(self):
        from eof.download import download_eofs

        # download all the misssed orbit files
        for record in self.df[self.df['orbitpath'].isna()].itertuples():
            #print (date, mission)
            orbitpath = download_eofs([record.datetime], [record.mission], save_dir=self.basedir)[0]
            #print ('orbitpath', orbitpath)
            self.df.loc[self.df.datetime == record.datetime,'orbitpath'] = orbitpath

    def set_dem(self, dem_filename):
        import os
        self.dem_filename = os.path.relpath(dem_filename,'.')
        return self

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small margin produces insufficient DEM not covers the defined area
    def download_dem(self, product='SRTM3', buffer_degrees=0.02, debug=False):
        import urllib.request
        import elevation
        import os
        import subprocess
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib

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
        # generate DEM for the area
        bounds = self.geoloc().dissolve().envelope.bounds.values[0]
        # show progress indicator
        def func():
            # left bottom right top
            elevation.clip(bounds=bounds,
                            product=product,
                            margin=str(buffer_degrees),
                            output=os.path.realpath(tif_filename))
            elevation.clean()
        with self.tqdm_joblib(notebook.tqdm(desc='Downloading', total=1)) as progress_bar:
            _ = joblib.Parallel(n_jobs=1)(joblib.delayed(func)() for i in [1])

        #if not os.path.exists(grd_filename):
        # convert to WGS84 ellipsoidal heights
        argv = ['gdalwarp', '-co', 'COMPRESS=DEFLATE',
                '-r', 'bilinear',
                '-s_srs', f'+proj=longlat +datum=WGS84 +no_defs +geoidgrids=./egm96_15.gtx',
                '-t_srs', '+proj=longlat +datum=WGS84 +no_def',
                '-overwrite',
                '-ot', 'Float32', '-of', 'NetCDF',
                'DEM_EGM96.tif', 'DEM_WGS84.nc']
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0:
            print ('download_dem', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('download_dem', stdout_data.decode('ascii'))

        self.dem_filename = grd_filename

    def backup(self, backup_dir):
        import os
        import shutil
    
        os.makedirs(backup_dir, exist_ok=True)
        for record in self.df.itertuples():
            for fname in [record.datapath, record.metapath, record.orbitpath]:
                shutil.copy2(fname, backup_dir)
        shutil.copy2(self.dem_filename, os.path.join(backup_dir, 'DEM_WGS84.nc'))
    
        return
    
    def set_master(self, master):
        if not master in self.df.index:
            raise Exception('Master image not found')
        self.master = master
        return self

    def get_master(self):
        return self.df.loc[[self.master]]

    def get_aligned(self, date=None):
        """
        Return selected aligned image or all the images (excluding master)
        """
        if date is None:
            idx = self.df.index.difference([self.master])
        else:
            idx = [date]
        return self.df.loc[idx]

#    def to_file(self, filename):
#        """
#        Save data.in file like to prep_data_linux.csh & prep_data.csh tools
#        """
#        if self.master is None:
#            raise Exception('Set master image first')
#        line = '\n'.join(self.get_master().apply(lambda row: f'{row.filename}:{row.orbitfile}', axis=1).values)
#        lines = '\n'.join(self.get_aligned().apply(lambda row: f'{row.filename}:{row.orbitfile}', axis=1).values)
#        with open(filename, 'wt') as f:
#            f.write(line+'\n'+lines+'\n')
#        return self

    def multistem_stem(self, dt=None):
        """
        Define stem and multistem using datetime    
        """
        from datetime import datetime

        # master datetime
        if dt is None:
            dt = self.df.loc[self.master, 'datetime']

        stem = f'S1_{dt.strftime("%Y%m%d_%H%M%S")}_F1'
        multistem = f'S1_{dt.strftime("%Y%m%d")}_ALL_F1'
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

#    def orbit(self, date=None):
#        """
#        Get orbit from XML scene annotation as single symbol 'A' | 'D'
#        """
#        
#        if date is None:
#            date = self.master
#       
#        filename = self.df.loc[date,'metapath']
#        doc = self.annotation(filename)
#        return doc['product']['generalAnnotation']['productInformation']['pass'][:1]

#    def geoloc(self, date=None):
#        """
#        Get GCPs from XML scene annotation as DataFrame
#        """
#       import pandas as pd
#        
#       if date is None:
#            date = self.master
#        
#        filename = self.df.loc[date,'metapath']
#        doc = self.annotation(filename)
#        #doc['geolocationGrid']
#        geoloc = doc['product']['geolocationGrid']['geolocationGridPointList']
#        # check data consistency
#        assert int(geoloc['@count']) == len(geoloc['geolocationGridPoint'])
#        geoloc_df = pd.DataFrame(geoloc['geolocationGridPoint']).applymap(lambda val : pd.to_numeric(val,errors='ignore'))
#        return geoloc_df

    def geoloc(self, filename=None):
        """
        Build approximate scene polygons using GCPs from XML scene annotation
        """
        from PRM import PRM
        import numpy as np
        import pandas as pd
        import geopandas as gpd
        import os

        if filename is None:
            filename = self.df.loc[self.master,'metapath']
        doc = self.annotation(filename)
        geoloc = doc['product']['geolocationGrid']['geolocationGridPointList']
        # check data consistency
        assert int(geoloc['@count']) == len(geoloc['geolocationGridPoint'])
        gcps = pd.DataFrame(geoloc['geolocationGridPoint']).applymap(lambda val : pd.to_numeric(val,errors='ignore'))
        # return approximate location as set of GCP
        return gpd.GeoDataFrame(gcps, geometry=gpd.points_from_xy(x=gcps.longitude, y=gcps.latitude))

#    # produce cropped frame using two pins
#    def geoloc_frame(self, pins):
#        """
#        Estimate framed area using two pins using Sentinel-1 GCPs with accuracy about 1 km.
#        The pins should be defined in any order as 1D or 2D array like to
#            [x1, y1, x2, y2] or [x2, y2, x1, y1] or [[x1, y1], [x2, y2]] or [[x2, y2], [x1, y1]]
#        The pins automatically reordered properly for ascending and descending orbits and returned in the true order.
#        """
#        import numpy as np
#        import geopandas as gpd
#        from shapely.geometry import LineString, Point, MultiPoint
#        from shapely.ops import split
#
#        def pin2line(pin, lons, lats):
#            coeffs = np.polyfit(lons, lats, 1)
#            poly = np.poly1d(coeffs)
#            dy = poly(pin.x) - pin.y
#            return LineString([Point(-180, poly(-180)-dy), Point( 180, poly( 180)-dy)])
#
#        date = self.master
#        
#        assert len(pins) == 4 or len(pins) == 2, 'Define two pins as two pairs of lat,lon coordinates'
#        orbit = self.orbit(self.master)
#        # convert to 1D array if needed
#        pins = np.array(pins).flatten()
#        assert len(pins) == 4, 'Define two pins as two pairs of lat,lon coordinates'
#        # swap pins if needed
#        if orbit == 'A' and pins[1] < pins[3]:
#            pin1 = Point(pins[:2])
#            pin2 = Point(pins[2:])
#        else:
#            pin1 = Point(pins[2:])
#            pin2 = Point(pins[:2])
#        df = self.geoloc()
#        area = MultiPoint(points=list(zip(df.lon, df.lat))).convex_hull
#        
#        line1 = pin2line(pin1,
#                         df[df.line==df.line.min()].lon,
#                         df[df.line==df.line.min()].lat)
#        line2 = pin2line(pin2,
#                         df[df.line==df.line.max()].lon,
#                         df[df.line==df.line.max()].lat)
#
#        diag = LineString([pin1,pin2])
#        geoms = [area.intersection(geom1).intersection(geom2) for geom1 in split(area, line1)
#                                                             for geom2 in split(area, line2)
#                if geom1.buffer(-1e-3).intersects(diag) and geom2.intersects(diag)
#                ]
#        assert len(geoms) > 0, 'Frame cannot be defined between the two pins. Hint: change the pins coordinates'
#        return gpd.GeoDataFrame({'name':['GCP','pin','pin','frame']},geometry=gpd.GeoSeries([area,pin1,pin2,geoms[0]]))

    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # small buffer produces incomplete area coverage and restricted NaNs
    # minimum buffer size: 8 arc seconds for 90 m DEM
    def get_dem(self, geoloc=False, buffer_degrees=0.02):
        import xarray as xr
        import os
        
        if self.dem_filename is None:
            raise Exception('Set DEM filename first')

        # open DEM file and find the elevation variable
        # because sometimes grid includes 'crs' or other variables
        dem = xr.open_dataset(self.dem_filename)
        assert 'lat' in dem.coords and 'lon' in dem.coords
        # define latlon array
        z_array_name = [data_var for data_var in dem.data_vars if len(dem.data_vars[data_var].coords)==2]
        assert len(z_array_name) == 1
        # extract the array and fill missed values by nan (mostly ocean area)
        dem = dem[z_array_name[0]].fillna(0)
        
        if geoloc is False:
            return dem
        
        bounds = self.geoloc().dissolve().envelope.bounds.values[0]
        #print ('xmin, xmax', xmin, xmax)
        return dem.sel(lat=slice(bounds[1]-buffer_degrees, bounds[3]+buffer_degrees),
                       lon=slice(bounds[0]-buffer_degrees, bounds[2]+buffer_degrees))

    def reframe(self, date, debug=False):
        """
        Estimate framed area using two pins using Sentinel-1 GCPs approximation.
        """
        from PRM import PRM
        import numpy as np
        import shapely
        import os

        df = self.df.loc[[date]]
        stem = self.multistem_stem(df['datetime'][0])[1]

        old_filename = os.path.join(self.basedir, f'{stem}')
        #print ('old_filename', old_filename)

        self.make_s1a_tops(date, debug)

        prm = PRM.from_file(old_filename+'.PRM')
        azi1 = prm.SAT_llt2rat([self.pins[0], self.pins[1], 0], precise=1)[1]
        azi2 = prm.SAT_llt2rat([self.pins[2], self.pins[3], 0], precise=1)[1]
        #print ('azi1', azi1, 'azi2', azi2)
        prm.shift_atime(azi1, inplace=True).update()

        # Working on bursts covering $azi1 ($ll1) - $azi2 ($ll2)...
        self.assemble_tops(date, azi1, azi2, debug)

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
            os.rename(old_filename+ext, new_filename+ext)

        # cleanup
        for fname in [old_filename+'.LED', old_filename+'.PRM']:
            if not os.path.exists(fname):
                continue
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

    def get_pins(self, pin_number=None):
        if pin_number is None:
            return self.pins
        elif pin_number == 1:
            return self.pins[:2]
        elif pin_number == 2:
            return self.pins[2:]

    def set_pins(self, pins):
        """
        Estimate framed area using two pins on Sentinel-1 GCPs approximation.
        The pins should be defined in any order like to
            [x1, y1, x2, y2]
            [[x1, y1], [x2, y2]]
            [[x1, y1], None]
            [None, [x2, y2]]
            [None, None]
        The pins automatically reordered properly for ascending and descending orbits and returned in the true order.
        """
        import numpy as np
        from shapely.geometry import Point

        error = False
        warning = False

        if pins is None:
            pins = [None, None]
        if len(pins) == 4:
            pins = np.array(pins).reshape(2,2)
        assert len(pins) == 2, 'Define two pins as two pairs of lat,lon coordinates where pin2 is upper pin1'
        #print ('pins', pins)
        pin1 = pins[0]
        pin2 = pins[1]
    
        df = self.df.loc[[self.master]]
        area = self.get_master()['geometry'].unary_union
        orbit = df['orbit'][0]

        # check the pins validity
        #geoloc = self.geoloc()
        llmin, ltmin, llmax, ltmax = self.geoloc().dissolve().envelope.bounds.values[0].round(3)
        #llmin, llmax = geoloc.longitude.min().round(3), geoloc.longitude.max().round(3)
        #ltmin, ltmax = geoloc.latitude.min().round(3), geoloc.latitude.max().round(3)
    
        if not np.all(pin1) is None and not area.intersects(Point(pin1[0], pin1[1])):
            print ('ERROR: pin1 lays outside of master frame. Move the pin or set it to None and try again.')
            error = True
        if not np.all(pin2) is None and not area.intersects(Point(pin2[0], pin2[1])):
            print ('ERROR: pin2 lays outside of master frame. Move the pin or set it to None and try again.')
            error = True
        
        # check pin1
        if np.all(pin1) is None and orbit == 'A':
            # use right bottom corner
            print ('NOTE: pin1 is not defined, master image corner coordinate will be used')
            warning = True
            pin1 = [llmax, ltmin]
        elif np.all(pin1) is None and orbit == 'D':
            # use right top corner
            print ('NOTE: pin1 is not defined, master image corner coordinate will be used')
            warning = True
            pin1 = [llmin, ltmin]

        # check pin2
        if np.all(pin2) is None and orbit == 'A':
            # use left top corner
            print ('NOTE: pin2 is not defined, master image corner coordinate will be used')
            warning = True
            pin2 = [llmin, ltmax]
        elif np.all(pin2) is None and orbit == 'D':
            # use left bottom corner
            print ('NOTE: pin2 is not defined, master image corner coordinate will be used')
            warning = True
            pin2 = [llmax, ltmax]

        # swap pins for Descending orbit
        if orbit == 'A':
            self.pins = [pin1[0], pin1[1], pin2[0], pin2[1]]
        else:
            self.pins = [pin2[0], pin2[1], pin1[0], pin1[1]]

        # pins are defined even is case of errors to have ability to plot them
        assert not error, 'ERROR: Please fix all the issues listed above to continue'

    def reframe_parallel(self, dates=None, n_jobs=-1, **kwargs):
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib
        import pandas as pd

        if dates is None:
            dates = self.df.index.unique().values

        # process all the scenes
        with self.tqdm_joblib(notebook.tqdm(desc='Reframing', total=len(dates))) as progress_bar:
            records = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.reframe)(date, **kwargs) for date in dates)

        self.df = pd.concat(records)

    def intf_ra2ll(self, grids):
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib
        import xarray as xr
        import numpy as np
        import os

        intf_ra2ll_file = os.path.join(self.basedir, 'intf_ra2ll.grd')

        matrix_ra2ll = xr.open_dataarray(intf_ra2ll_file)

        def ra2ll(grid):
            return xr.DataArray(np.where(matrix_ra2ll>=0, grid.values.reshape(-1)[matrix_ra2ll], np.nan),
                coords=matrix_ra2ll.coords)

        # process single 2D raster
        if len(grids.dims) == 2:
            return ra2ll(grids)

        # process a set of 2D rasters
        with self.tqdm_joblib(notebook.tqdm(desc='Geocoding', total=len(grids))) as progress_bar:
            grids_ll = joblib.Parallel(n_jobs=-1)(joblib.delayed(ra2ll)(grids[item]) for item in range(len(grids)))
        grids_ll = xr.concat(grids_ll, dim=grids.dims[0])
    
        # add coordinates from original grids
        for coord in grids.coords:
            if coord in ['y', 'x']:
                continue
            grids_ll[coord] = grids[coord]

        return grids_ll

    def intf_ra2ll_matrix(self, intf_grids):
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np
        import os

        trans_dat_file = os.path.join(self.basedir, 'trans.dat')
        trans_ra2ll_file = os.path.join(self.basedir, 'trans_ra2ll.grd')
        intf_ra2ll_file = os.path.join(self.basedir, 'intf_ra2ll.grd')
        
        # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        trans = np.fromfile(trans_dat_file, dtype=np.float64).reshape([-1,5])
        lon_min, lon_max = trans[:,3].min(),trans[:,3].max()
        lat_min, lat_max = trans[:,4].min(),trans[:,4].max()

        trans_ra2ll = xr.open_dataarray(trans_ra2ll_file)

        intf_ys, intf_xs = xr.broadcast(intf_grids[0].y, intf_grids[0].x)
        intf_yxs = np.stack([intf_ys.values.reshape(-1),intf_xs.values.reshape(-1)], axis=1)
        trans_yxs = np.stack([trans[:,1],trans[:,0]], axis=1)

        tree = cKDTree(intf_yxs, compact_nodes=False, balanced_tree=False)
        distance_limit = np.max([intf_grids[0].y.diff('y')[0], intf_grids[0].x.diff('x')[0]])
        d, inds = tree.query(trans_yxs, k = 1, distance_upper_bound=distance_limit, workers=8)

        # single integer index mask
        intf2trans = np.where(~np.isinf(d), inds, -1)
        # produce the same output array
        intf_ra2ll = xr.zeros_like(trans_ra2ll).rename('intf_ra2ll')
        intf_ra2ll.values = np.where(trans_ra2ll>=0, intf2trans[trans_ra2ll], -1)
        # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
        intf_ra2ll.attrs['node_offset'] = 1
        # save to NetCDF file
        intf_ra2ll.to_netcdf(intf_ra2ll_file, encoding={'intf_ra2ll': self.compression})

    def ra2ll(self):
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np
        import os

        trans_dat_file = os.path.join(self.basedir, 'trans.dat')
        trans_ra2ll_file = os.path.join(self.basedir, 'trans_ra2ll.grd')
        
        if os.path.exists(trans_ra2ll_file):
            os.remove(trans_ra2ll_file)

        # trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        trans = np.fromfile(trans_dat_file, dtype=np.float64).reshape([-1,5])
        lon_min, lon_max = trans[:,3].min(),trans[:,3].max()
        lat_min, lat_max = trans[:,4].min(),trans[:,4].max()

        #dem = xr.open_dataset(in_dem_gridfile)
        #dem = self.get_dem(geoloc=True)
        dem = self.get_dem(geoloc=True).sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        
        trans_latlons = np.stack([trans[:,4],trans[:,3]], axis=1)
        dem_lats, dem_lons = xr.broadcast(dem.lat,dem.lon)
        dem_latlons = np.stack([dem_lats.values.reshape(-1),dem_lons.values.reshape(-1)], axis=1)

        tree = cKDTree(trans_latlons, compact_nodes=False, balanced_tree=False)
        distance_limit = np.max([dem.lat.diff('lat')[0],dem.lon.diff('lon')[0]])
        d, inds = tree.query(dem_latlons, k = 1, distance_upper_bound=distance_limit, workers=8)

        # produce the same output array as dataset to be able to add global attributes
        trans_ra2ll = xr.zeros_like(dem).rename('trans_ra2ll')
        trans_ra2ll.values = np.where(~np.isinf(d), inds, -1).reshape(dem.shape)
        # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
        #trans_ra2ll.attrs['node_offset'] = 1
        # save to NetCDF file
        trans_ra2ll.to_netcdf(trans_ra2ll_file, encoding={'trans_ra2ll': self.compression})

    def topo_ra(self, method='cubic'):
        import xarray as xr
        import os

        trans_dat_file = os.path.join(self.basedir, 'trans.dat')
        topo_ra_file = os.path.join(self.basedir, 'topo_ra.grd')

        # cleanup before createing the new files
        for filename in [trans_dat_file, topo_ra_file]:
            if os.path.exists(filename):
                os.remove(filename)
        
        self.PRM().topo_ra(self.get_dem(geoloc=True),
                           trans_dat_tofile=trans_dat_file,
                           topo_ra_tofile=topo_ra_file,
                           method=method)
        
        # build DEM coordinates transform matrix
        self.ra2ll()

    # replacement for gmt grdfilter ../topo/dem.grd -D2 -Fg2 -I12s -Gflt.grd
    def get_topo_llt(self, degrees, geoloc=True):
        import xarray as xr
        import numpy as np

        # add buffer around the cropped area for borders interpolation
        dem_area = self.get_dem(geoloc=geoloc)
        ny = int(np.round(degrees/dem_area.lat.diff('lat')[0]))
        nx = int(np.round(degrees/dem_area.lon.diff('lon')[0]))
        #print ('DEM decimation','ny', ny, 'nx', nx)
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

    def assemble_tops(self, date, azi_1, azi_2, debug=False):
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

        df = self.df.loc[[date]]
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
        stem = self.multistem_stem(df['datetime'][0])[1]

        # round values and convert to strings
        azi_1 = np.round(azi_1).astype(int).astype(str)
        azi_2 = np.round(azi_2).astype(int).astype(str)

        argv = ['assemble_tops', azi_1, azi_2] + datapaths + [stem]
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('assemble_tops', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('assemble_tops', stdout_data.decode('ascii'))

        return

    def ext_orb_s1a(self, stem, date=None, debug=False):
        import os
        import subprocess

        if date is None or date == self.master:
            df = self.get_master()
        else:
            df = self.df.loc[[date]]

        orbit = os.path.relpath(df['orbitpath'][0], self.basedir)
    
        argv = ['ext_orb_s1a', f'{stem}.PRM', orbit, stem]
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('ext_orb_s1a', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('ext_orb_s1a', stdout_data.decode('ascii'))

        return
    
    # produce LED and PRM in basedir
    # when date=None work on master image
    def make_s1a_tops(self, date=None, mode=0, rshift_fromfile=None, ashift_fromfile=None, debug=False):
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
            date = self.master
            # for master image mode should be 1
            mode = 1
        df = self.df.loc[[date]]
    
        xmlfile = os.path.relpath(df['metapath'][0], self.basedir)
        datafile = os.path.relpath(df['datapath'][0], self.basedir)
        stem = self.multistem_stem(df['datetime'][0])[1]
    
        argv = ['make_s1a_tops', xmlfile, datafile, stem, str(mode)]
        if rshift_fromfile is not None:
            argv.append(rshift_fromfile)
        if ashift_fromfile is not None:
            argv.append(ashift_fromfile)
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('make_s1a_tops', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('make_s1a_tops', stdout_data.decode('ascii'))
    
        self.ext_orb_s1a(stem, date, debug=debug)
    
        return

    # aligning for master image
    def stack_ref(self, debug=False):
        import xarray as xr
        import numpy as np
        import os
        from PRM import PRM

        err, warn = self.validate()
        #print ('err, warn', err, warn)
        assert not err and not warn, 'ERROR: Please fix all the issues listed above to continue'

        master_line = list(self.get_master().itertuples())[0]
        #print (master_line)

        # for master image
        multistem, stem = self.multistem_stem(master_line.datetime)
        path_stem = os.path.join(self.basedir, stem)
        path_multistem = os.path.join(self.basedir, multistem)

        # generate PRM, LED, SLC
        self.make_s1a_tops(debug=debug)

        PRM.from_file(path_stem + '.PRM')\
            .set(input_file = path_multistem + '.raw')\
            .update(path_multistem + '.PRM', safe=True)

        self.ext_orb_s1a(multistem, debug=debug)

        # recalculate after ext_orb_s1a
        earth_radius = PRM.from_file(path_multistem + '.PRM')\
            .calc_dop_orb(inplace=True).update().get('earth_radius')

    # aligning for secondary image
    def stack_rep(self, date=None, degrees=12.0/3600, debug=False):
        import xarray as xr
        import numpy as np
        import os
        from PRM import PRM

        err, warn = self.validate()
        #print ('err, warn', err, warn)
        assert not err and not warn, 'ERROR: Please fix all the issues listed above to continue'

        master_line = list(self.get_master().itertuples())[0]
        multistem, stem = self.multistem_stem(master_line.datetime)
        #print (master_line)

        # define master image parameters
        master = self.PRM().sel('earth_radius').set(stem=stem, multistem=multistem)

        # prepare coarse DEM for alignment
        # 12 arc seconds resolution is enough, for SRTM 90m decimation is 4x4
        topo_llt = self.get_topo_llt(degrees=degrees)
        #topo_llt.shape

        line = list(self.get_aligned(date).itertuples())[0]
        multistem, stem = self.multistem_stem(line.datetime)
        #print (line)

        # define relative filenames for PRM
        stem_prm    = os.path.join(self.basedir, stem + '.PRM')
        mstem_prm   = os.path.join(self.basedir, multistem + '.PRM')
        master_prm  = os.path.join(self.basedir, master.get("stem") + '.PRM')
        mmaster_prm = os.path.join(self.basedir, master.get("multistem") + '.PRM')

        # TODO: define 1st image for line, in the example we have no more
        tmp_da = 0

        # generate PRM, LED
        self.make_s1a_tops(date, debug=debug)

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
        r_grd.to_netcdf(stem_prm[:-4]+'_r.grd')

        a_grd = self.offset2shift(a_xyz, rmax, amax)
        a_grd.to_netcdf(stem_prm[:-4]+'_a.grd')

        # generate the image with point-by-point shifts
        # note: it removes calc_dop_orb parameters from PRM file
        # generate PRM, LED
        self.make_s1a_tops(date=line.Index, mode=1,
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

        self.ext_orb_s1a(multistem, date=line.Index, debug=debug)

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

    def stack_parallel(self, dates=None, n_jobs=-1, **kwargs):
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib
    
        if dates is None:
            dates = list(self.get_aligned().index)

        # prepare master image
        self.stack_ref()

        # prepare secondary images
        with self.tqdm_joblib(notebook.tqdm(desc='Aligning', total=len(dates))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.stack_rep)(date, **kwargs) for date in dates)

    def intf(self, pair, **kwargs):
        from PRM import PRM
        import os

        # extract dates from pair
        date1, date2 = pair

        prm_ref = self.PRM(date1)
        prm_rep = self.PRM(date2)

        #print ('SBAS intf kwargs', kwargs)
        prm_ref.intf(prm_rep, basedir=self.basedir, **kwargs)

    def intf_parallel(self, pairs, n_jobs=-1, **kwargs):
        import pandas as pd
        import numpy as np
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib
        from joblib.externals import loky
        import os

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        # this way does not work properly for long interferogram series
        #with self.tqdm_joblib(notebook.tqdm(desc='Interferograms', total=len(pairs))) as progress_bar:
        #    joblib.Parallel(n_jobs=-1)(joblib.delayed(self.intf)(pair, **kwargs) for pair in pairs)

        # start a set of jobs together but not more than available cpu cores at once
        if n_jobs == -1:
            n_jobs = joblib.cpu_count()
        n_chunks = int(np.ceil(len(pairs)/n_jobs))
        chunks = np.array_split(pairs, n_chunks)
        with notebook.tqdm(desc='Interferograms', total=len(pairs)) as pbar:
            for chunk in chunks:
                loky.get_reusable_executor(kill_workers=True).shutdown(wait=True)
                with joblib.parallel_backend('loky', n_jobs=len(chunk), inner_max_num_threads=1):
                    joblib.Parallel()(joblib.delayed(self.intf)(pair, **kwargs) for pair in chunk)
                pbar.update(len(chunk))

        # build radar coordinates transformation matrix
        self.intf_ra2ll_matrix(self.open_grids(pairs, 'phasefilt'))
        
        # build stack mask
        filename_mask = os.path.join(self.basedir, 'mask.grd')
        mask = self.open_grids(pairs, 'mask').mean('pair')
        mask.to_netcdf(filename_mask, encoding={'z': self.compression})
        
        # build stack coherence
        filename_corr = os.path.join(self.basedir, 'corr.grd')
        corr = self.open_grids(pairs, 'corr').mean('pair')
        corr.to_netcdf(filename_mask, encoding={'z': self.compression})
        

    def baseline_table(self):
        import pandas as pd
        
        prm_ref = self.PRM()
        data = []
        for date in self.df.index:
            prm_rep = self.PRM(date)
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

    @staticmethod
    def baseline_triplets_analysis(baseline_triplets, grids):
        import pandas as pd
        import numpy as np

        data = []
        for triplet in baseline_triplets.itertuples():
            #print (triplet)
            pair_a = f'{triplet.A} {triplet.B}'
            pair_b = f'{triplet.B} {triplet.C}'
            pair_c = f'{triplet.A} {triplet.C}'
            #print (pair_a, pair_b, pair_c)
            grid_a = grids.sel(pair=pair_a)
            grid_b = grids.sel(pair=pair_b)
            grid_c = grids.sel(pair=pair_c)
            grid_abc = (grid_a + grid_b - grid_c)
            diff_abc = float(np.round(grid_abc.mean() % np.pi,3))
            delta_pi = int(np.round(grid_abc.mean() / np.pi,0))

            a = grid_a.values.ravel()
            b = grid_b.values.ravel()
            c = grid_c.values.ravel()
            mask = np.isnan(a)|np.isnan(b)|np.isnan(c)

            corr_a = (np.corrcoef(a[~mask], b[~mask])[1,0].round(2))
            corr_b = (np.corrcoef(b[~mask], c[~mask])[1,0].round(2))
            corr_c = (np.corrcoef(a[~mask], c[~mask])[1,0].round(2))

            data.append({'A': triplet.A, 'B': triplet.B, 'C': triplet.C,
                         'corr_AB': corr_a, 'corr_BC': corr_b, 'corr_CA' :corr_c, 'delta_pi': delta_pi})

        return (pd.DataFrame(data).sort_values(['A', 'B', 'C']))

    def PRM(self, date=None, multi=True):
        from PRM import PRM
        import os

        if date is None or date == self.master:
            line = self.get_master()
        else:
            line = self.get_aligned(date)
        #print (line)
        multistem, stem = self.multistem_stem(line.datetime[0])
        if multi:
            stem = multistem
        filename = os.path.join(self.basedir, f'{stem}.PRM')
        #print (filename)
        return PRM.from_file(filename)

    def unwrap_parallel(self, pairs, n_jobs=-1, **kwargs):
        import pandas as pd
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        with self.tqdm_joblib(notebook.tqdm(desc='Unwrapping', total=len(pairs))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.unwrap)(pair, **kwargs) for pair in pairs)

    # -s for SMOOTH mode and -d for DEFO mode when DEFOMAX_CYCLE should be defined in the configuration
    # DEFO mode (-d) and DEFOMAX_CYCLE=0 is equal to SMOOTH mode (-s)
    def unwrap(self, pair, threshold, conf=None, func=None, debug=False):
        import xarray as xr
        import numpy as np
        import os
        import subprocess

        if conf is None:
            conf = self.PRM().snaphu_config(defomax=0)

        # extract dates from pair
        date1, date2 = pair

        basename = os.path.join(self.basedir, f'{date1}_{date2}_').replace('-','')
        #print ('basename', basename)

        # invalid pixels mask
        mask = xr.open_dataarray(basename + 'mask.grd')
    
        phase = xr.open_dataarray(basename + 'phasefilt.grd')
        phase_in = basename + 'unwrap.phase'
        phase.where(~np.isnan(phase),0).values.tofile(phase_in)

        corr = xr.open_dataarray(basename + 'corr.grd')
        # mask invalid pixels on correlation matrix
        corr = (corr*mask).where(corr>=threshold)
        corr_in = basename + 'unwrap.corr'
        corr.where(~np.isnan(corr),0).values.tofile(corr_in)
        
        # TEST
        if os.path.exists(basename + 'corr.tmp.grd'):
            os.remove(basename + 'corr.tmp.grd')
        corr.where(~np.isnan(corr),0).to_netcdf(basename + 'corr.tmp.grd', encoding={'z': self.compression})

        unwrap_out = basename + 'unwrap.out'

        argv = ['snaphu', phase_in, str(phase.shape[1]), '-c', corr_in,
                '-f', '/dev/stdin', '-o', unwrap_out, '-d']
        if debug:
            argv.append('-v')
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='ascii', bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=conf)
        if len(stderr_data) > 0 and debug:
            print ('snaphu', stderr_data)
        if debug:
            print ('snaphu', stdout_data)

        # read results
        values = np.fromfile(unwrap_out, dtype=np.float32).reshape(phase.shape)
        #values = np.frombuffer(stdout_data, dtype=np.float32).reshape(mask.shape)
        # mask invalid pixels on unwrapped results
        unwrap = mask * xr.DataArray(values, phase.coords, name='z')
        # apply user-defined function for post-processing
        if func is not None:
            unwrap = func(corr, unwrap)
        if os.path.exists(basename + 'unwrap.grd'):
            os.remove(basename + 'unwrap.grd')
        # mask again when user-defined function applied
        (mask * unwrap).to_netcdf(basename + 'unwrap.grd', encoding={'z': self.compression})

        for tmp_file in [phase_in, corr_in, unwrap_out]:
            #print ('tmp_file', tmp_file)
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

    #gmt grdmath unwrap_mask.grd $wavel MUL -79.58 MUL = los.grd
    def los_displacement_mm(self, unwraps):
        # constant is negative to make LOS = -1 * range change
        # constant is (1000 mm) / (4 * pi)
        scale = -79.58 * self.PRM().get('radar_wavelength')
        return scale*unwraps

    # returns all grids in basedir by mask or grids by dates and name
    def open_grids(self, pairs, name, geocode=False, mask=False, func=None, crop_valid=False):
        import pandas as pd
        import xarray as xr
        import numpy as np
        import os

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values
        else:
            pairs = np.asarray(pairs)

        #if pairs is None:
        #    # stack by filepath for xr.open_mfdataset
        #    def preprocess_dirname(ds):
        #        pair = os.path.basename(ds.encoding['source'])[:17]
        #        #print (ds.encoding['source'], '->', pair)
        #        return ds.assign(pair=pair)
        #    filenames = os.path.join(self.basedir, f'*_{name}.grd')
        #    ds = xr.open_mfdataset(filenames, concat_dim='pair', combine='nested',
        #                             preprocess=preprocess_dirname)['z']
        #    return ds

        das = []
        if len(pairs.shape) == 1:
            for date in sorted(pairs):
                filename = os.path.join(self.basedir, f'{name}_{date}.grd'.replace('-',''))
                #print (date, filename)
                da = xr.open_dataarray(filename)
                if func is not None:
                    da = func(da)
                if mask:
                    da = da*self.open_grid('mask')
                das.append(da.expand_dims('date'))
            das = xr.concat(das, dim='date')
            das['date'] = sorted(pairs)
        elif len(pairs.shape) == 2:
            for pair in pairs:
                filename = os.path.join(self.basedir, f'{pair[0]}_{pair[1]}_{name}.grd'.replace('-',''))
                #print (filename)
                da = xr.open_dataarray(filename)
                if func is not None:
                    da = func(da)
                if mask:
                    da = da*self.open_grid('mask')
                das.append(da.expand_dims('pair'))
            das = xr.concat(das, dim='pair')
            das['pair'] = [f'{pair[0]} {pair[1]}' for pair in pairs]
            das['ref']  = xr.DataArray([pair[0] for pair in pairs], dims='pair')
            das['rep']  = xr.DataArray([pair[1] for pair in pairs], dims='pair')
        else:
            raise Exception('Use single or two columns Pandas dataset or array as "pairs" argument')

        if geocode:
            das = self.intf_ra2ll(das)

        # crop NaNs
        if crop_valid:
            dims = [dim for dim in das.dims if dim != 'pair']
            assert len(dims) == 2, 'ERROR: interferogram should be 2D array'
            da = das.min('pair')
            indexer = {}
            for dim in dims:
                da = da.dropna(dim=dim, how='all')
                dim_min, dim_max = da[dim].min().item(), da[dim].max().item()
                indexer[dim] = slice(dim_min, dim_max)
            #print ('indexer', indexer)
            return das.loc[indexer]

        return das

    def open_grid(self, name, geocode=False, mask=False, func=None):
        import pandas as pd
        import xarray as xr
        import os

        filename = os.path.join(self.basedir, f'{name}.grd')
        da = xr.open_dataarray(filename)
        if func is not None:
            da = func(da)
        if mask:
            #da = da*self.open_grid('mask')
            da = da * self.open_grid('mask').interp_like(da, method='nearest')
        if geocode:
            return self.intf_ra2ll(da)
        return da

    #intf.tab format:   unwrap.grd  corr.grd  ref_id  rep_id  B_perp 
    def intftab(self, baseline_pairs):
        import numpy as np
        import datetime

        outs = []
        for line in baseline_pairs.itertuples():
            #print (line)
            ref = line.ref_date.replace('-','')
            jref = datetime.datetime.strptime(line.ref_date, '%Y-%m-%d').strftime('%Y%j')
            rep = line.rep_date.replace('-','')
            jrep = datetime.datetime.strptime(line.rep_date, '%Y-%m-%d').strftime('%Y%j')
            bperp = np.round(line.rep_baseline - line.ref_baseline, 2)
            outs.append(f'{ref}_{rep}_unwrap.grd {ref}_{rep}_corr.grd {jref} {jrep} {bperp}')
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

    def sbas(self, baseline_pairs, smooth=0, debug=False):
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
        dem = self.get_dem()
        geoloc = self.geoloc()
        prm = self.PRM()

        #N=$(wc -l intf.in   | cut -d ' ' -f1)
        #S=$(wc -l scene.tab | cut -d ' ' -f1)

        N = len(baseline_pairs)
        S = len(np.unique(list(baseline_pairs['ref_date']) + list(baseline_pairs['rep_date'])))

        #bounds = self.geoloc().dissolve().envelope.bounds.values[0]
        lon0 = geoloc.longitude.mean()
        lat0 = geoloc.latitude.mean()
        elevation0 = float(dem.sel(lat=lat0, lon=lon0, method='nearest'))
        #print ('coords',lon0, lat0, elevation0)
        _,_,_,look_E,look_N,look_U = prm.SAT_look([lon0, lat0, elevation0])[0]
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
                str(N), str(S), str(xdim), str(ydim), '-smooth', str(smooth),
                '-wavelength', str(wavelength), '-incidence', str(incidence), '-range', str(rng), '-rms', '-dem']
        if debug:
            print (' '.join(argv))
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=self.basedir, encoding='ascii')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('sbas', stderr_data)
        if len(stdout_data) > 0 and debug:
            print ('sbas', stdout_data)

        # fix output grid filenames
        for date in np.unique(np.concatenate([baseline_pairs.ref_date,baseline_pairs.rep_date])):
            jdate = datetime.datetime.strptime(date, '%Y-%m-%d').strftime('%Y%j')
            date = date.replace('-','')
            filename1 = os.path.join(self.basedir, f'disp_{jdate}.grd')
            filename2 = os.path.join(self.basedir, f'disp_{date}.grd')
            if os.path.exists(filename1):
                os.rename(filename1, filename2)
            #print (jdate, date)

        return

#filelist = SBAS('raw_orig').set_master(MASTER)
#filelist.df
#filelist.get_master()
#filelist.get_aligned()
#filelist.to_file('raw/data.in.new')
