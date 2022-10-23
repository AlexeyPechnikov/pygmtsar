#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to scan *xml files and orbits and make data.in file like to prep_data_linux.csh & prep_data.csh tools
from .tqdm_joblib import tqdm_joblib
from .gmtsar import gmtsar
from .snaphu import snaphu
from .orbits import orbits
from .dem import dem
from .landmask import landmask
from .geocode import geocode

class SBAS(tqdm_joblib, gmtsar, snaphu, orbits, dem, landmask, geocode):

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

        # see https://github.com/mobigroup/gmtsar/issues/8
        df = df.sort_values(by=['date', 'subswath']).set_index('date')\
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
        import numpy as np
        import shapely
        import os
        from pygmtsar import PRM

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
        #if len(subswaths) == 1:
        #    # build geo transform matrices for interferograms
        #    self.transforms(subswaths[0], pairs)

    def merge(self, pair, grid, debug=False):
        import os
        from pygmtsar import PRM

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
        #self.topo_ra_parallel()
        #self.transforms(pairs)

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
        import os
        from pygmtsar import PRM

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

    # use the function for open_grids() and save_grids()
    def get_filenames(self, subswath, pairs, name, add_subswath=True):
        import pandas as pd
        import numpy as np
        import os

        # define subswath when subswath=None and check otherwise
        subswath = self.get_subswath(subswath)
    
        if pairs is None:
            # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
            pass
        elif isinstance(pairs, pd.DataFrame):
            pairs = pairs.values
        else:
            pairs = np.asarray(pairs)

        if add_subswath == True:
            prefix = f'F{subswath}_'
        else:
            prefix = ''

        filenames = []
        if pairs is None:
            # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
            filename = os.path.join(self.basedir, f'{prefix}{name}.grd')
            return filename
        elif len(pairs.shape) == 1:
            # read all the grids from files
            for date in sorted(pairs):
                filename = os.path.join(self.basedir, f'{prefix}{name}_{date}.grd'.replace('-',''))
                filenames.append(filename)
        elif len(pairs.shape) == 2:
            # read all the grids from files
            for pair in pairs:
                filename = os.path.join(self.basedir, f'{prefix}{pair[0]}_{pair[1]}_{name}.grd'.replace('-',''))
                filenames.append(filename)
        return filenames

    def save_grids(self, grids, name, func=None, add_subswath=True, n_jobs=1, **kwargs):
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        import os

        subswaths = self.get_subswaths()
        # pipeline before merging subswaths is well defined and we do not need to use this function
        assert len(subswaths) == 1, 'ERROR: use SBAS.merge() to merge multiple subswaths before'
    
        assert len(grids.dims) in [2, 3], 'ERROR: supported 2D and 3D arrays only'
        if len(grids.dims) ==3 :
            grids = [da for da in grids]
    
        # special case for a single grid
        if not isinstance(grids, list):
            filename = self.get_filenames(None, None, name, add_subswath=add_subswath)
            grids.astype(np.float32).rename(name)\
                .to_netcdf(filename, encoding={name: self.compression}, engine=self.engine)
            return
    
        def preprocess(subswath, da):
            if func is not None:
                if isinstance(func, list):
                    for f in func:
                        da = f(da)
                else:
                    da = func(da)
            # check 3rd dimension
            if 'date' in da[0].coords:
                # "date" for SBAS results for an example
                pair = da.date.item()
            elif 'pair' in grids[0].coords:
                # "pair" for interferograms for an example
                pair = da.pair.item().split(' ')
            # the function returns a set of filenames
            filename = self.get_filenames(subswath, [pair], name, add_subswath=add_subswath)[0]
            #print ('filename', filename)
            if os.path.exists(filename):
                os.remove(filename)
            da.astype(np.float32).rename(name)\
                .to_netcdf(filename, encoding={name: self.compression}, engine=self.engine)
            return

        # test
        #for subswath in subswaths:
        #    for grid in grids:
        #        preprocess(subswath, grid)
        #return

        # process all the grids
        description = 'Saving' if func is None else 'Processing and Saving'
        with self.tqdm_joblib(tqdm(desc=description, total=len(grids)*len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(preprocess)(subswath, grid) \
                                                     for subswath in subswaths for grid in grids)

    # returns all grids in basedir by mask or grids by dates and name
    # Backward-compatible open_grids() returns list of grids fot the name or a single grid for a single subswath
    def open_grids(self, pairs, name, geocode=False, inverse_geocode=False,  mask=None, func=None,
                   crop_valid=False, add_subswath=True, chunks=None, n_jobs=-1):
        import pandas as pd
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm
        import joblib

        assert not(geocode and inverse_geocode), 'ERROR: Only single geocoding option can be applied'

        if chunks is None:
            chunks = self.chunksize

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
            filenames = self.get_filenames(subswath, pairs=pairs, name=name, add_subswath=add_subswath)
        
            das = []
            if pairs is None:
                # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
                #print ('filename', filename)
                da = xr.open_dataarray(filenames, engine=self.engine, chunks=chunks)
                das  = postprocess(da, subswath)
            elif len(pairs.shape) == 1:
                # read all the grids from files
                for filename in filenames:
                    #print (date, filename)
                    da = xr.open_dataarray(filename, engine=self.engine, chunks=chunks)
                    das.append(da)

                # post-processing on a set of 2D rasters
                with self.tqdm_joblib(tqdm(desc='Loading', total=len(das))) as progress_bar:
                    das = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(postprocess)(da, subswath) for da in das)

                # prepare to stacking
                das = [da.expand_dims('date') for da in das]

                # allow stack to be extended to largest 1st dimension size
                # to be sure all code work well for this case
                # so user is able to load grids by his own way
                das = xr.concat(das, dim='date')
                das['date'] = sorted(pairs)
            elif len(pairs.shape) == 2:
                # read all the grids from files
                for filename in filenames:
                    #print (filename)
                    da = xr.open_dataarray(filename, engine=self.engine, chunks=chunks)
                    das.append(da)

                # post-processing on a set of 2D rasters
                with self.tqdm_joblib(tqdm(desc='Loading', total=len(das))) as progress_bar:
                    das = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(postprocess)(da, subswath) for da in das)

                # prepare to stacking
                das = [da.expand_dims('pair') for da in das]

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

    def detrend(self, dataarray, fit_intercept=True, fit_dem=True, fit_coords=True,
                resolution_meters=90, debug=False):
        """
        Detrend unwrapped interferogram in radar coordinates, see for details
            https://github.com/gmtsar/gmtsar/issues/98
            https://github.com/gmtsar/gmtsar/issues/411
        """
        import xarray as xr
        import numpy as np
        import dask
        from sklearn.linear_model import LinearRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
    
        def postprocessing(out):
            return out.astype(np.float32).rename('detrend')

        # check the simplest case
        assert fit_intercept or fit_dem or fit_coords, 'All the detrending options disable, function does nothing'

        # check simple case
        if fit_intercept and not fit_dem and not fit_coords:
            if debug:
                print ('DEBUG: Remove mean value only')
            return postprocessing(dataarray - dataarray.mean())

        # input grid can be too large
        decimator = self.pixel_decimator(resolution_meters=resolution_meters, grid=dataarray, debug=debug)
        # decimate
        dataarray_dec = decimator(dataarray)

        # topography grid required to fit_dem option only
        if fit_dem:
            if debug:
                print ('DEBUG: Interpolate topography on the data grid')
            topo_ra = self.get_topo_ra()
            #topo = topo.reindex_like(unwraps[0], method='nearest')
            # use xr.zeros_like to prevent the target grid coordinates modifying
            topo = topo_ra.reindex_like(dataarray, method='nearest')
            # check chunks
            if debug:
                print ('DEBUG: regrid to resolution in meters', resolution_meters)
            # decimate
            topo_dec  = decimator(topo)
        else:
            topo = topo_dec = None

        # lazy calculations are useless below
        def data2fit(data, grid):
            y = data.values.reshape(-1)
            nanmask = np.isnan(y)
            # prepare regression variable
            Y = y[~nanmask]

            if fit_coords or fit_dem:
                # prepare coordinates for X regression variable
                yy, xx = xr.broadcast(data.y, data.x)
                ys = yy.values.reshape(-1)[~nanmask]
                xs = xx.values.reshape(-1)[~nanmask]

            if fit_dem:
                # prepare topography for X regression variable
                zs = grid.values.reshape(-1)[~nanmask]
                zys = zs*ys
                zxs = zs*xs

            if fit_dem and fit_coords:
                if debug:
                    print ('DEBUG: Detrend topography and coordinates')
                X = np.column_stack([zys, zxs, ys, xs, zs])
            elif fit_dem:
                if debug:
                    print ('DEBUG: Detrend topography only')
                X = np.column_stack([zys, zxs, zs])
            elif fit_coords:
                if debug:
                    print ('DEBUG: Detrend coordinates only')
                X = np.column_stack([ys, xs])
            return Y, X, nanmask

        # build prediction model with or without plane removal (fit_intercept)
        regr = make_pipeline(StandardScaler(), LinearRegression(fit_intercept=fit_intercept))
        Y, X, _ = data2fit(dataarray_dec, topo_dec)
        regr.fit(X, Y)
    
        # TODO: calculate for chunks
        Y, X, nanmask = data2fit(dataarray, topo)
        model = xr.full_like(dataarray, np.nan).compute()
        if debug:
            print ('DEBUG: model', model)
        model.data.reshape(-1)[~nanmask] = regr.predict(X)
        return postprocessing(dataarray - model)

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
        netcdf_compression = self.compression.copy()
        netcdf_compression['chunksizes'] = (1, self.chunksize, self.chunksize)

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
                         engine=self.engine,
                         compute=True)
            return chunk_filename
    
        # process all the chunks
        with self.tqdm_joblib(tqdm(desc='Computing', total=ys*xs)) as progress_bar:
            filenames = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(func)(iy, ix) \
                                                     for iy in range(ys) for ix in range(xs))
    
        # rebuild the datasets to user-friendly format
        das = [xr.open_dataarray(f, engine=self.engine, chunks=self.chunksize) for f in filenames]
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
                        encoding={'displacement': self.compression},
                        engine=self.engine)

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
