# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------

class S1():

    @staticmethod
    def scan_slc(datadir, orbit=None, mission=None, subswath=None, polarization=None):
        """
        Initialize an instance of the Stack class.

        Parameters
        ----------
        datadir : str
            The directory containing the data files.
        dem_filename : str, optional
            The filename of the DEM (Digital Elevation Model) WGS84 NetCDF file. Default is None.
        landmask_filename : str, optional
            The filename of the landmask WGS84 NetCDF file. Default is None.
        orbit : str, optional
            Filter for orbit direction. Use 'A' for Ascending, 'D' for Descending, or None for no filter. Default is None.
        mission : str, optional
            Filter for mission name. Use 'S1A' for Sentinel-1A, 'S1B' for Sentinel-1B, or None for no filter. Default is None.
        subswath : int, optional
            Filter for subswath number. Use a single or sequential numbers 1, 2, 3, 12, 23, 123, or None for no filter. Default is None.
        polarization : str, optional
            Filter for polarization. Use 'VV', 'VH', 'HH', 'HV', or None for no filter. Default is None.

        Examples
        --------
        Initialize an Stack object with the data directory 'data':
        stack = S1('data')
        """
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

        #text2date('V20171110T225942'), text2date('20171117t145927')
        def text2date(text, as_date=True):
            """
            Convert a text string in the format 'VYYYYMMDDTHHMMSS' or 'YYYYMMDDTHHMMSS' to a datetime object or date object.

            Parameters
            ----------
            text : str
                The text string to convert.
            as_date : bool, optional
                If True, return a date object. If False, return a datetime object. Default is True.

            Returns
            -------
            datetime.date or datetime.datetime
                The converted date or datetime object.
            """
            from datetime import datetime
        
            date_fmt = '%Y%m%dT%H%M%S'
            date_str = text.replace('V','').replace('t','T')
            dt = datetime.strptime(date_str, date_fmt)
            if as_date:
                return dt.date()
            return dt

        def pattern2paths(pattern):
            path_pattern = os.path.join(datadir, '**', pattern)
            paths = glob(path_pattern, recursive=True)
            return paths

        assert orbit is None or orbit=='A' or orbit=='D', \
            'ERROR: use symbol A (Ascending) or D (Descending) for orbit filter'
        assert mission is None or mission=='S1A' or mission=='S1B', \
            'ERROR: use name S1A or S1B for mission filter'
        assert subswath is None or subswath in [1,2,3,12,23,123], \
            'ERROR: use a single or sequential numbers 1, 2, 3, 12, 23, 123 for subswath filter'
        assert polarization is None or polarization in ['VV','VH','HH','HV'], \
            'ERROR: use VV or VH or HH or HV for polarization filter'

        if polarization is None:
            polarization = '??'
        if subswath is None:
            subswath  = '?'
        else:
            subswath = f'[{subswath}]'
        # filter mission
        if mission is not None:
            path_pattern = f'{mission.lower()}-iw{subswath}-slc-{polarization.lower()}-*'
        else:
            path_pattern = f's1?-iw{subswath}-slc-{polarization.lower()}-*'
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
        dts = [text2date(name.split('-')[4],False) for name in datanames]
        #print ('filedatetimes', dts)

        ds = [dt.date() for dt in dts]
        #print ('filedates', ds)

        df = pd.DataFrame({'date':[str(d) for d in ds], 'datetime': dts, 'datapath': datapaths, 'metapath': metapaths})
        #print ('df', df)
        assert len(df), f'Scenes not found'

        # filter mission and always ignore approximate RESORB orbits to download precise POEORB when possible
        if mission is not None:
            orbit_path_pattern = f'{mission.upper()}_OPER_AUX_???ORB_OPOD_*.EOF'
        else:
            orbit_path_pattern = 'S1?_OPER_AUX_???ORB_OPOD_*.EOF'
        orbitpaths = pattern2paths(orbit_path_pattern)
        #print ('orbitpaths', orbitpaths)
        orbitnames = [os.path.splitext(os.path.split(path)[-1])[0] for path in orbitpaths]
        if orbitpaths:
            orbit_dates = [(text2date(name.split('_')[-2]), text2date(name.split('_')[-1])) for name in orbitnames]
            orbits = dict(zip(orbit_dates, orbitpaths))
            #print ('orbits', orbits)
            # look for as precise (from date-1 day to date+1 day) as restituted orbits (from date to date or date-1 to date)
            orbits = [orbits.get((date-oneday, date+oneday)) or
                      orbits.get((date-oneday, date)) or
                      orbits.get((date, date)) for date in ds]
            #print ('fileorbits', fileorbits)
            df['orbitpath'] = orbits
        else:
            df['orbitpath'] = None

        # add some calculated properties
        df['subswath'] = [int(filename.split('-')[1][-1]) for filename in datanames]
        df['mission'] = [filename.split('-')[0].upper() for filename in datanames]
        df['polarization'] = [filename.split('-')[3].upper() for filename in datanames]
    
        # read approximate locations
        #geolocs = [shapely.geometry.MultiPoint(self.geoloc(path).geometry).minimum_rotated_rectangle for path in metapaths]
        #print ('geolocs', geolocs)
        #df = gpd.GeoDataFrame(df, geometry=geolocs)
        def geoloc2bursts(metapath):
            """
            Read approximate bursts locations
            """
            from shapely.geometry import LineString, Polygon, MultiPolygon
            df = S1.get_geoloc(S1.read_annotation(metapath))
            # this code line works for a single scene
            #lines = df.groupby('line')['geometry'].apply(lambda x: LineString(x.tolist()))
            # more complex code is required for stitched scenes processing with repeating 'line' series
            df['line_change'] = df['line'].diff().ne(0).cumsum()
            # single-point lines possible for stitched scenes
            grouped_lines = df.groupby('line_change')['geometry'].apply(lambda x: LineString(x.tolist()) if len(x) > 1 else None)
            lines = grouped_lines.reset_index(drop=True)
            #bursts = [Polygon([*line1.coords, *line2.coords[::-1]]) for line1, line2 in zip(lines[:-1], lines[1:])]
            # to ignore None for single-point lines
            bursts = []
            prev_line = None
            for line in lines:
                if line is not None and prev_line is not None:
                    bursts.append(Polygon([*prev_line.coords, *line.coords[::-1]]))
                prev_line = line
            return MultiPolygon(bursts)
        bursts = [geoloc2bursts(path) for path in metapaths]
        df = gpd.GeoDataFrame(df, geometry=bursts)

        # define orbit directions
        orbits = [S1.read_annotation(path)['product']['generalAnnotation']['productInformation']['pass'][:1] for path in metapaths]
        df['orbit'] = orbits
        # filter orbits
        if orbit is not None:
            df = df[df.orbit == orbit]
        #print ('df', df)
        assert len(df), f'Scenes not found for the defined orbit {orbit}'

        # see https://github.com/mobigroup/gmtsar/issues/8
        df = df.sort_values(by=['date', 'subswath']).set_index('date')\
            [['datetime','orbit','mission','polarization','subswath','datapath','metapath','orbitpath','geometry']]

        # Validate the DataFrame to check for any issues or inconsistencies.
        # we can't merge together scenes from different missions
        missions = df.groupby('date')['mission'].unique().values
        missions = [len(mission) for mission in missions if len(mission)>1]
        if not len(missions) == 0:
            raise ValueError('ERROR: Found multiple scenes for a single date from different missions')
        # note: df.unique() returns unsorted values so it would be 21 instead of expected 12
        subswaths = int(''.join(map(str,np.unique(df.subswath))))
        if not int(subswaths) in [1, 2, 3, 12, 23, 123]:
            raise ValueError(f'ERROR: Subswhats list {subswaths} incorrect. Allowed a single or sequential subswath numbers 1, 2, 3, 12, 23, 123')
        if not len(df.orbit.unique()) <= 1:
            raise ValueError('ERROR: Only single orbit processing supported. Use any one ascending or descending')
        if not len(df.index.unique()) >= 2:
            raise ValueError('ERROR: Two or more scenes required')
        daily_scenes = df.groupby(['date', 'subswath'])['datetime'].count().values.max()
        if daily_scenes > 1:
            print ('NOTE: Found multiple scenes for a single day, use function Stack.reframe() to stitch the scenes')

        return df

    @staticmethod
    def read_annotation(filename):
        """
        Return the XML scene annotation as a dictionary.

        Parameters
        ----------
        filename : str
            The filename of the XML scene annotation.

        Returns
        -------
        dict
            The XML scene annotation as a dictionary.
        """
        import xmltodict

        with open(filename) as fd:
            # fix wrong XML tags to process cropped scenes
            # GMTSAR assemble_tops.c produces malformed xml
            # https://github.com/gmtsar/gmtsar/issues/354
            doc = xmltodict.parse(fd.read().replace('/></','></'))
        return doc

    @staticmethod
    def get_geoloc(annotation):
        """
        Build approximate scene polygons using Ground Control Points (GCPs) from XML scene annotation.

        Parameters
        ----------
        filename : str, optional
            The filename of the XML scene annotation. If None, print a note and return an empty DataFrame. Default is None.

        Returns
        -------
        geopandas.GeoDataFrame
            A GeoDataFrame containing the approximate scene polygons.

        annotation = S1.read_annotation(filename)
        S1.get_geoloc(annotation)
        """
        import numpy as np
        import pandas as pd
        import geopandas as gpd
        import os
    
        geoloc = annotation['product']['geolocationGrid']['geolocationGridPointList']
        # check data consistency
        assert int(geoloc['@count']) == len(geoloc['geolocationGridPoint'])
        # Google Colab wrapper
        if float(pd.__version__[:3]) > 2.0:
            gcps = pd.DataFrame(geoloc['geolocationGridPoint']).map(lambda val : pd.to_numeric(val,errors='ignore'))
        else:
            gcps = pd.DataFrame(geoloc['geolocationGridPoint']).applymap(lambda val : pd.to_numeric(val,errors='ignore'))
        # return approximate location as set of GCP
        return gpd.GeoDataFrame(gcps, geometry=gpd.points_from_xy(x=gcps.longitude, y=gcps.latitude))
