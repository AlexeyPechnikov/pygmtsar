#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# Wrapper to scan *xml files and orbits and make data.in file like to prep_data_linux.csh & prep_data.csh tools
from .SBAS_incidence import SBAS_incidence
from .PRM import PRM

class SBAS(SBAS_incidence):

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

    def PRM(self, subswath=None, date=None, multi=True, singleswath=False):
        """
        multi=True/False - open multistem or stem file
        singleswath=False/True - open a single-digit subswath instead of a multi-digit (merged) one
            single-digit subswath exists always while multi-digit exists only for interferogram pair references
        """
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

#    def make_gaussian_filter(self, range_dec, azi_dec, wavelength, debug=False):
#        """
#        Wrapper for PRM.make_gaussian_filter() and sonamed command line tool. Added for development purposes only.
#        """
#        import numpy as np
#
#        gauss_dec, gauss_string = self.PRM().make_gaussian_filter(range_dec, azi_dec, wavelength, debug=debug)
#        coeffs = [item for line in gauss_string.split('\n') for item in line.split('\t') if item != '']
#        # x,y dims order
#        shape = np.array(coeffs[0].split(' ')).astype(int)
#        # y,x dims order
#        matrix = np.array(coeffs[1:]).astype(float).reshape((shape[1],shape[0]))
#        return (gauss_dec, matrix)
