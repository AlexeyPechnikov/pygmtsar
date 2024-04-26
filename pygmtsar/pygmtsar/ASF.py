# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .tqdm_joblib import tqdm_joblib
from .S1 import S1

class ASF(tqdm_joblib):
    import pandas as pd
    from datetime import timedelta

    template_url = 'https://datapool.asf.alaska.edu/SLC/S{satellite}/{scene}.zip'
    template_safe = '*.SAFE/*/{mission}-iw{subswath}-slc-{polarization}-*'
    # check for downloaded scene files
    template_scene = 'S1?_IW_SLC__1S??_{start}_*.SAFE/*/s1?-iw{subswath_lowercase}-slc-{polarization_lowercase}-{start_lowercase}-*'
    # URL of the HTML page
    # see https://asf.alaska.edu/data-sets/sar-data-sets/sentinel-1/sentinel-1-data-and-imagery/
    # https://s1qc.asf.alaska.edu
    poeorb_url = 'https://s1qc.asf.alaska.edu/aux_poeorb/'
    resorb_url = 'https://s1qc.asf.alaska.edu/aux_resorb/'
    # resorb - 69108
    # poeorb - 9811
    #pattern_orbit  = r'S1\w_OPER_AUX_\w{3}ORB_OPOD_\d{8}T\d{6}_V\d{8}T\d{6}_\d{8}T\d{6}.EOF'
    pattern_orbit  = r'S1\w_OPER_AUX_(POE|RES)ORB_OPOD_\d{8}T\d{6}_V\d{8}T\d{6}_\d{8}T\d{6}.EOF'
    # see _select_orbit.py in sentineleof package
    #Orbital period of Sentinel-1 in seconds
    #T_ORBIT = (12 * 86400.0) / 175.0
    offset_start = timedelta(seconds=(12 * 86400.0) // 175.0 + 60)
    offset_end = timedelta(seconds=300)

    def __init__(self, username=None, password=None):
        import asf_search
        import getpass
        if username is None:
            username = getpass.getpass('Please enter your ASF username and press Enter key:')
        if password is None:
            password = getpass.getpass('Please enter your ASF password and press Enter key:')
        self.username = username
        self.password = password

    def _get_asf_session(self):
        import asf_search
        return asf_search.ASFSession().auth_with_creds(self.username, self.password)

    def download_scenes(self, basedir, scenes, subswaths, polarization='VV', session=None, n_jobs=4, skip_exist=True):
        import pandas as pd
        import numpy as np
        import asf_search
        import fnmatch
        import joblib
        from tqdm.auto import tqdm
        import os
        import re
        import glob
        from datetime import datetime, timedelta
        import warnings
        # supress asf_search 'UserWarning: File already exists, skipping download'
        warnings.filterwarnings("ignore", category=UserWarning)

        # create the directory if needed
        os.makedirs(basedir, exist_ok=True)

        # collect all the existing files once
        files = os.listdir(basedir)
        #print ('files', len(files))

        # skip existing scenes
        if skip_exist:
            # collect all the existing files once
            files = glob.glob('**', root_dir=basedir, recursive=True)
            #print ('files', len(files))
            # check scenes folders only
            #scenes_missed = np.unique([scene for scene in scenes if f'{scene}.SAFE' not in files])
            scenes_missed = []
            for scene in scenes:
                scene_path = '/'.join([scene + '.SAFE'] + self.template_scene.split('/')[1:])
                #print ('scene_path', scene_path)
                for subswath in str(subswaths):
                    pattern = scene_path.format(
                                         subswath_lowercase = subswath.lower(),
                                         polarization_lowercase = polarization.lower(),
                                         start_lowercase = '*')
                    #print ('pattern', pattern)
                    matching = [filename for filename in files if fnmatch.fnmatch(filename, pattern)]
                    exts = [os.path.splitext(fname)[1] for fname in matching]
                    if '.tiff' in exts and '.xml'in exts:
                        pass
                    else:
                        scenes_missed.append(scene)
    
        else:
            # process all the defined scenes
            scenes_missed = scenes
        #print ('scenes_missed', len(scenes_missed))

        # do not use internet connection, work offline when all the scenes and orbits already available
        if len(scenes_missed) == 0:
            return

        def get_url(scene):
            return self.template_url.format(satellite=scene[2:3], scene=scene)

        def get_patterns(subswaths, polarization, mission):
            #print (f'get_patterns: {subswaths}, {polarization}, {mission}')
            assert len(mission) == 3 and mission[:2]=='S1', \
                f'ERROR: mission name is invalid: {mission}. Expected names like "S1A", "S1B", etc.'
            outs = []
            for subswath in str(subswaths):
                pattern = self.template_safe.format(mission=mission.lower(),
                                                    subswath=subswath,
                                                    polarization=polarization.lower())
                outs.append(pattern)
            return outs

        def download_scene(scene, subswaths, polarization, basedir, session):
            # define scene zip url
            url = get_url(scene)
            #print (f'download_scene: {url}, {patterns}, {basedir}, {session}')
            patterns = get_patterns(subswaths, polarization, mission=scene[:3])
            #print ('patterns', patterns)
            with asf_search.remotezip(url, session) as remotezip:
                filenames = remotezip.namelist()
                #print ('filenames', filenames)
                for pattern in patterns:
                    #print ('pattern', pattern)
                    matching = [filename for filename in filenames if fnmatch.fnmatch(filename, pattern)]
                    for filename in matching:
                        filesize = remotezip.getinfo(filename).file_size
                        #print (filename, filesize)
                        fullname = os.path.join(basedir, filename)
                        if os.path.exists(fullname) and os.path.getsize(fullname) == filesize:
                            #print (f'pass {fullname}')
                            pass
                        else:
                            #print (f'download {fullname}')
                            # create the directory if needed
                            try:
                                os.makedirs(os.path.dirname(fullname), exist_ok=True)
                                if os.path.exists(fullname + '.tmp'):
                                    os.remove(fullname + '.tmp')
                                with open(fullname + '.tmp', 'wb') as file:
                                    file.write(remotezip.read(filename))
                                assert os.path.getsize(fullname + '.tmp') == filesize, \
                                    f'ERROR: Downloaded incomplete scene content'
                                os.rename(fullname + '.tmp', fullname)
                            except Exception as e:
                                print(e)
                                raise
                            finally:
                                if os.path.exists(fullname + '.tmp'):
                                    os.remove(fullname + '.tmp')
                            #remotezip.extract(filename, basedir)

        # prepare authorized connection
        if session is None:
            session = self._get_asf_session()

        # download scenes
        with self.tqdm_joblib(tqdm(desc='ASF Downloading Sentinel-1 SLC:', total=len(scenes_missed))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(download_scene)\
                                    (scene, subswaths, polarization, basedir, session) for scene in scenes_missed)

        # parse processed scenes and convert to dataframe
        #print ('scenes', len(scenes))
        scenes_downloaded = pd.DataFrame(scenes_missed, columns=['scene'])
        # return the results in a user-friendly dataframe
        #scenes_downloaded['scene'] = scenes_downloaded.scene\
        #                             .apply(lambda name: self.template_url.format(satellite=name[2:3], scene=name))
        return scenes_downloaded

    def download(self, *args, **kwarg):
        print ('NOTE: Function is deprecated. Use ASF.download_scenes() and ASF.download_orbits().')

    # https://asf.alaska.edu/datasets/data-sets/derived-data-sets/sentinel-1-bursts/
    def download_bursts(self, basedir, bursts, session=None, n_jobs=4, skip_exist=True):
        import rioxarray as rio
        import xmltodict
        import pandas as pd
        import asf_search
        import joblib
        from tqdm.auto import tqdm
        import os
        import glob
        from datetime import datetime, timedelta
        import warnings
        # supress asf_search 'UserWarning: File already exists, skipping download'
        warnings.filterwarnings("ignore", category=UserWarning)
    
        # create the directory if needed
        os.makedirs(basedir, exist_ok=True)
    
        # skip existing bursts
        if skip_exist:
            bursts_missed = []
            for burst in bursts:
                #print (burst)
                burst_parts = burst.split('_')
                subswath = burst_parts[2]
                polarization = burst_parts[4]
                start = burst_parts[3]
                #print ('start', start, 'subswath', subswath, 'polarization', polarization)
                template = self.template_scene.format(start=start,
                                                     subswath_lowercase = subswath.lower(),
                                                     polarization_lowercase = polarization.lower(),
                                                     start_lowercase = start.lower())
                #print ('template', template)
                files = glob.glob(template, root_dir=basedir)
                exts =[ os.path.splitext(file)[-1] for file in files]
                if sorted(exts) == ['.tiff', '.xml']:
                    #print ('pass')
                    pass
                else:
                    bursts_missed.append(burst)
        else:
            # process all the defined scenes
            bursts_missed = bursts
        #print ('bursts_missed', len(bursts_missed))
        # do not use internet connection, work offline when all the scenes and orbits already available
        if len(bursts_missed) == 0:
            return
    
        def download_burst(result, basedir, session):
            properties = result.geojson()['properties']
            #print (properties)
            polarization = properties['polarization']
            #print ('polarization', polarization)
            subswath = properties['burst']['subswath']
            #print ('subswath', subswath)
            # fake image number unique per polarizations
            # the real image number is sequential between all the polarizations
            # this trick allows to detect scene files without manifest downloading
            imageNumber = '00' + subswath[-1:]
            #print ('imageNumber', imageNumber)
            # define new scene name
            scene = properties['url'].split('/')[3]
             #print ('scene', scene)
            # use fake start and stop times to follow the burst naming
            start_time = properties['sceneName'].split('_')[3]
            start_time_dt = datetime.strptime(start_time, '%Y%m%dT%H%M%S')
            stop_time_dt = (start_time_dt  + timedelta(seconds=3))
            stop_time = stop_time_dt.strftime('%Y%m%dT%H%M%S')
            #print ('start_time, stop_time', start_time, stop_time)
            scene_parts = scene.split('_')
            #scene_parts[5] = startTime.replace('-','').replace(':','')[:len(scene_parts[5])]
            scene_parts[5] = start_time
            #scene_parts[6] = stopTime.replace('-','').replace(':','')[:len(scene_parts[6])]
            scene_parts[6] = stop_time
            scene = '_'.join(scene_parts)
            #print ('scene', scene)
            scene_dir = os.path.join(basedir, scene + '.SAFE')
            # define burst name
            burst = '-'.join(['s1a-iw1-slc-vv'] + scene_parts[5:-1] + [imageNumber]).lower()
            # 's1a-iw1-slc-vv-20240314t130744-20240314t130747-052978-0669c6-005'
            #print ('burst', burst)
    
            # create the directories if needed
            tif_dir = os.path.join(scene_dir, 'measurement')
            xml_dir = os.path.join(scene_dir, 'annotation')
            # save annotation using the burst and scene names
            xml_file = os.path.join(xml_dir, burst + '.xml')
            #rint ('xml_file', xml_file)
            tif_file = os.path.join(tif_dir, burst + '.tiff')
            #print ('tif_file', tif_file)
            for dirname in [scene_dir, tif_dir, xml_dir]:
                os.makedirs(dirname, exist_ok=True)
    
            # download tif
            # properties['bytes'] is not an accurate file size but it looks about 40 kB smaller
            if os.path.exists(tif_file) and os.path.getsize(tif_file) >= properties['bytes']:
                #print (f'pass {tif_file}')
                pass
            else:
                #print ('YYY', os.path.getsize(tif_file), properties['bytes'])
                # remove potentially incomplete file if needed
                if os.path.exists(tif_file):
                    os.remove(tif_file)
                # download burst tif file and save using the burst and scene names
                #result.download(os.path.dirname(tif_file), filename=os.path.basename(tif_file))
                result.download(scene_dir, filename=os.path.basename(tif_file))
                # check if we can open the downloaded file without errors
                tmp_file = os.path.join(scene_dir, os.path.basename(tif_file))
                with rio.open_rasterio(tmp_file) as raster:
                    raster.load()
                    os.rename(tmp_file, tif_file)
        
            # download xml
            if os.path.exists(xml_file) and os.path.getsize(xml_file) > 0:
                #print (f'pass {xml_file}')
                pass
            else:
                # get the file name
                basename = os.path.basename(properties['additionalUrls'][0])
                #print ('basename', '=>', basename)
                manifest_file = os.path.join(scene_dir, basename)
                # remove potentially incomplete manifest file if needed
                if os.path.exists(manifest_file):
                    os.remove(manifest_file)
                # download and process manifest file even when it exists but is not processed to annotation xml
                asf_search.download_urls(urls=properties['additionalUrls'], path=scene_dir, session=session)
                # parse xml
                with open(manifest_file, 'r') as file:
                    xml_content = file.read()
                subswathidx = int(subswath[-1:]) - 1 
                annotation = xmltodict.parse(xml_content)['burst']['metadata']['product'][subswathidx]['content']
                #imageNumber = annotation['adsHeader']['imageNumber']
                #print ('imageNumber', imageNumber)
                # filter GCPs for only the selected burst properties['burst']['burstIndex']
                # filter GCPs by time
                geoloc = annotation['geolocationGrid']['geolocationGridPointList']
                # check data consistency
                assert int(geoloc['@count']) == len(geoloc['geolocationGridPoint'])
                # to produce 2 lines instead of single use time offset
                startTime = properties['startTime']
                #print ('startTime', startTime)
                stopTime  = properties['stopTime']
                #print ('stopTime', stopTime)
                startTime_dt_buffer = datetime.strptime(startTime, '%Y-%m-%dT%H:%M:%S.%fZ') + timedelta(seconds=-3)
                startTime_buffer = startTime_dt_buffer.strftime('%Y-%m-%dT%H:%M:%S.%f') + 'Z'
                stopTime_dt_buffer = datetime.strptime(stopTime, '%Y-%m-%dT%H:%M:%S.%fZ') + timedelta(seconds=0)
                stopTime_buffer = stopTime_dt_buffer.strftime('%Y-%m-%dT%H:%M:%S.%f') + 'Z'
                #print ('stopTime_buffer', stopTime_buffer)
                gcps = []
                for gcp in geoloc['geolocationGridPoint']:
                    if gcp['azimuthTime'] >= startTime_buffer and gcp['azimuthTime'] <= stopTime_buffer:
                        #print (gcp)
                        gcps.append(gcp)
                geoloc['@count'] = str(len(gcps))
                geoloc['geolocationGridPoint'] = gcps
                annotation['geolocationGrid']['geolocationGridPointList'] = geoloc
                with open(xml_file, 'w') as file:
                    file.write(xmltodict.unparse({'product': annotation}, pretty=True, indent='  '))
                # remove processed manifest file
                if os.path.exists(manifest_file):
                    os.remove(manifest_file)
    
        # prepare authorized connection
        if session is None:
            session = self._get_asf_session()
    
        with tqdm(desc=f'ASF Downloading Bursts Catalog', total=1) as pbar:
            results = asf_search.granule_search(bursts_missed)
            pbar.update(1)
    
        # download bursts
        with self.tqdm_joblib(tqdm(desc='ASF Downloading Sentinel-1 Bursts', total=len(bursts_missed))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(download_burst)\
                                    (result, basedir, session) for result in results)
    
        # parse processed bursts and convert to dataframe
        bursts_downloaded = pd.DataFrame(bursts_missed, columns=['burst'])
        # return the results in a user-friendly dataframe
        return bursts_downloaded
