# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .S1_base import S1_base

class S1_slc(S1_base):

    template_burst = '*_*_?W?/{type}/S1_*_?W?_*_??_????-BURST'
    template_orbit = 'S1?_OPER_AUX_???ORB_OPOD_*.EOF'

    @staticmethod
    def scan(datadir):
        """
        Scans the specified directory for Sentinel-1 SLC (Single Look Complex) data and filters it based on the provided parameters.
    
        Parameters
        ----------
        datadir : str
            The directory containing the data files.
        
        Returns
        -------
        pandas.DataFrame
            A DataFrame containing metadata about the found burst, including their paths and other relevant properties.
    
        Raises
        ------
        ValueError
            If the bursts contain inconsistencies, such as mismatched .tiff and .xml files, or if invalid filter parameters are provided.
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
            path_pattern = os.path.join(datadir, pattern)
            #print ('path_pattern', path_pattern)
            paths = glob(path_pattern, recursive=True)
            return paths

        datapaths = pattern2paths(S1_slc.template_burst.format(type='measurement') + '.tiff')
        #print ('datapaths', datapaths)
        paths = [os.path.basename(os.path.dirname(os.path.dirname(p))) for p in datapaths]
        metapaths = pattern2paths(S1_slc.template_burst.format(type='annotation') + '.xml')
        #print ('metapaths', metapaths)
        noisepaths = pattern2paths(S1_slc.template_burst.format(type='noise') + '.xml')
        #print ('noisepaths', noisepaths)
        calibpaths = pattern2paths(S1_slc.template_burst.format(type='calibration') + '.xml')
        #print ('calibpaths', calibpaths)

        datanames = [os.path.splitext(os.path.basename(path))[0] for path in datapaths]
        #print ('datanames', datanames)
        metanames = [os.path.splitext(os.path.basename(path))[0] for path in metapaths]
        #print ('metanames', metanames)
        noisenames = [os.path.splitext(os.path.basename(path))[0] for path in noisepaths]
        #print ('noisenames', noisenames)
        calibnames = [os.path.splitext(os.path.basename(path))[0] for path in calibpaths]
        #print ('calibnames', calibnames)

        datas = dict(zip(datanames, datapaths))
        metas = dict(zip(metanames, metapaths))
        noises = dict(zip(noisenames, noisepaths))
        calibs = dict(zip(calibnames, calibpaths))

        # define the same order when and only when the names are the same
        assert sorted(datanames) == sorted(metanames), 'Found inconsistent set of .tiff and .xml files'
        # reorder paths using the same order
        datapaths = [datas[name] for name in datanames]
        metapaths = [metas[name] for name in metanames]
        assert sorted(datanames) == sorted(noisenames), 'Found inconsistent set of .xml files'
        assert sorted(noisenames) == sorted(calibnames), 'Found inconsistent set of calibration .xml files'
        noisepaths = [noises[name] for name in noisenames]
        calibpaths = [calibs[name] for name in calibnames]

        # points to datadir and extensions tiff, xml
        #print ('filenames', filenames)
        dts = [text2date(name.split('_')[3],False) for name in datanames]
        #print ('filedatetimes', dts)

        ds = [dt.date() for dt in dts]
        #print ('filedates', ds)

        df = pd.DataFrame({
                            'date':[str(d) for d in ds],
                            'datetime': dts,
                            'burst': datanames,
                            'fullBurstID': paths
                           })
        assert len(df), f'Scenes not found'

        bursts = [S1_slc.geoloc2bursts(path) for path in metapaths]
        df = gpd.GeoDataFrame(df, geometry=bursts)

        # define orbit directions
        flightDirections = [S1_slc.read_annotation(path)['product']['generalAnnotation']['productInformation']['pass'].upper() for path in metapaths]
        df['flightDirection'] = flightDirections

        modes = [S1_slc.read_annotation(path)['product']['adsHeader']['mode'] for path in metapaths]
        df['beamModeType'] = modes
        
        missions = [S1_slc.read_annotation(path)['product']['adsHeader']['missionId'] for path in metapaths]
        df['mission'] = missions
        
        subswaths = [S1_slc.read_annotation(path)['product']['adsHeader']['swath'] for path in metapaths]
        df['subswath'] = subswaths
        
        polarizations = [S1_slc.read_annotation(path)['product']['adsHeader']['polarisation'] for path in metapaths]
        df['polarization'] = polarizations
        
        # ((absolute_orbit - 73) % 175) + 1
        absolute_orbits = [S1_slc.read_annotation(path)['product']['adsHeader']['absoluteOrbitNumber'] for path in metapaths]
        df['pathNumber'] = [((int(absolute_orbit) - 73) % 175) + 1 for absolute_orbit in absolute_orbits]

        #return df

        # always ignore approximate RESORB orbits to download precise POEORB when possible
        orbitpaths = pattern2paths(S1_slc.template_orbit)
        #print ('orbitpaths', orbitpaths)
        orbitnames = [os.path.split(path)[-1] for path in orbitpaths]
        #print ('orbitnames', orbitnames)
        orbit_dates = [
            (text2date(start), text2date(end))
            for name in orbitnames
            for _, start, end in [os.path.splitext(name)[0].rsplit('_', 2)]
        ]
        orbits = dict(zip(orbit_dates, orbitnames))
        #print ('orbits', orbits)
        # look for as precise (from date-1 day to date+1 day) as restituted orbits (from date to date or date-1 to date)
        orbits = [orbits.get((date-oneday, date+oneday)) or
                  orbits.get((date-oneday, date)) or
                  orbits.get((date, date)) for date in ds]
        #print ('fileorbits', fileorbits)
        df['orbit'] = orbits

        # see https://github.com/mobigroup/gmtsar/issues/8
        df = df.sort_values(by=['fullBurstID', 'date']).set_index(['fullBurstID', 'date'])\
            [['datetime','mission','beamModeType','flightDirection','pathNumber','polarization','subswath','burst', 'orbit','geometry']]

        # specify the scan path
        df['path'] = datadir
        return df

    @staticmethod
    def geoloc2bursts(metapath):
        """
        Read approximate bursts locations
        """
        from shapely.geometry import LineString, Polygon, MultiPolygon
        df = S1_slc.get_geoloc(S1_slc.read_annotation(metapath))
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
    
        gcps = pd.DataFrame(geoloc['geolocationGridPoint'])
        # convert to numeric values excluding azimuthTime & slantRangeTime
        for column in gcps.columns[2:]:
            gcps[column] = pd.to_numeric(gcps[column])

        # return approximate location as set of GCP
        return gpd.GeoDataFrame(gcps, geometry=gpd.points_from_xy(x=gcps.longitude, y=gcps.latitude))
