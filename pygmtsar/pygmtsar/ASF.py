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

class ASF(tqdm_joblib):
    from datetime import timedelta

    template_url = 'https://datapool.asf.alaska.edu/SLC/S{satellite}/{scene}.zip'
    template_safe = '*.SAFE/*/s1{mission}-iw{subswath}-slc-{polarization}-*'
    # URL of the HTML page
    # see https://asf.alaska.edu/data-sets/sar-data-sets/sentinel-1/sentinel-1-data-and-imagery/
    # https://s1qc.asf.alaska.edu
    poeorb_url = 'https://s1qc.asf.alaska.edu/aux_poeorb/'
    resorb_url = 'https://s1qc.asf.alaska.edu/aux_resorb/'
    # resorb - 69108
    # poeorb - 9811
    #pattern_orbit  = r'S1\w_OPER_AUX_\w{3}ORB_OPOD_\d{8}T\d{6}_V\d{8}T\d{6}_\d{8}T\d{6}.EOF'
    pattern_orbit  = r'S1\w_OPER_AUX_(POE|RES)ORB_OPOD_\d{8}T\d{6}_V\d{8}T\d{6}_\d{8}T\d{6}.EOF'
    template_orbit = '{mission}_OPER_AUX_(POE|RES)ORB_OPOD_\d+T\d+_V{date_start}T\d+_{date_end}T\d+.EOF'
    # see _select_orbit.py in sentineleof package
    #Orbital period of Sentinel-1 in seconds
    #T_ORBIT = (12 * 86400.0) / 175.0
    orbit_start_offset = timedelta(seconds=(12 * 86400.0) // 175.0 + 60)
    orbit_end_offset = timedelta(seconds=300)
    
    def __init__(self, username=None, password=None):
        import asf_search
        import getpass
        if username is None:
            username = getpass.getpass('Please enter your ASF username and press Enter key:')
        if password is None:
            password = getpass.getpass('Please enter your ASF password and press Enter key:')
        self.username = username
        self.password = password

    def _get_orbits(self, scenes, url, session):
        import pandas as pd
        import re
        from datetime import datetime

        #print ('url', url)
        # get orbits HTML list
        responce = session.get(url)
        if responce.status_code == 200:
            lines = responce.text.splitlines()
            # TODO: cleanup code
            orbits = [re.search(self.pattern_orbit, line).group(0) if re.search(self.pattern_orbit, line) else None for line in lines]
            orbits = pd.DataFrame(orbits, columns=['orbit']).dropna()
        else:
            print('ERROR: Failed to fetch the web page. Status code: {response.status_code}')
            # return the input as is
            return scenes

        orbits['mission']    = orbits['orbit'].apply(lambda name: name[:3])
        orbits['time']       = orbits['orbit'].apply(lambda name: datetime.strptime(name.split('_')[5],      '%Y%m%dT%H%M%S'))
        orbits['time_start'] = orbits['orbit'].apply(lambda name: datetime.strptime(name.split('_')[6][1:],  '%Y%m%dT%H%M%S'))
        orbits['time_end']   = orbits['orbit'].apply(lambda name: datetime.strptime(name[:-4].split('_')[7], '%Y%m%dT%H%M%S'))

        for ind, scene in enumerate(scenes.itertuples()):
            #print (scene)
            if scene.orbit is not None:
                continue
            # look for the orbit file
            orbit = orbits[(orbits.mission == scene.mission)&\
                           (orbits.time_start <= scene.time_start - self.orbit_start_offset)&\
                           (orbits.time_end   >= scene.time_end   + self.orbit_end_offset  )].sort_values('time')
            # add the recent orbit
            orbit = orbit.tail(1).orbit.item() if len(orbit) >= 1 else None
            if orbit is None:
                continue
            scenes.at[ind, 'orbit'] = url + orbit
        return scenes

    """
    find . -type f -name '*.tiff' -exec basename {} .tiff \; \
        | awk -F'-' -v OFS='_' '{print toupper($1), "IW_SLC__1SDV", toupper($5), toupper($6), $7, toupper($8), $9"X.SAFE"}' \
        | xargs -I {} mkdir -p {}
    """
    def download(self, basedir, scenes, subswaths, polarization='VV', mission='?', n_jobs=8, skip_exist=True):
        import pandas as pd
        import numpy as np
        import asf_search
        import fnmatch
        import joblib
        from tqdm.auto import tqdm
        import os
        import re
        from datetime import datetime, timedelta
        import warnings
        # supress asf_search 'UserWarning: File already exists, skipping download'
        warnings.filterwarnings("ignore", category=UserWarning)

        # create the directory if needed
        os.makedirs(basedir, exist_ok=True)

        # check if scenes and orbits are already downloaded, no internet connection required
        files = os.listdir(basedir)
        scenes_exist = []
        for scene in scenes:
            # check scene
            if not f'{scene}.SAFE' in files:
                continue
            # check orbit
            mission = scene[:3]
            dt = datetime.strptime(scene.split('_')[5], '%Y%m%dT%H%M%S')
            day        = dt.strftime('%Y%m%d')
            day_before = (dt + timedelta(days=-1)).strftime('%Y%m%d')
            day_after  = (dt + timedelta(days=1)).strftime('%Y%m%d')
            pattern    = self.template_orbit.format(mission=mission, date_start=day_before, date_end=day_after)
            filenames  = sorted(filter(re.compile(pattern).match, files))
            if len(filenames) == 0:
                pattern   = self.template_orbit.format(mission=mission, date_start=day, date_end=day)
                filenames = sorted(filter(re.compile(pattern).match, files))
            if len(filenames) == 0:
                continue
            # scene and orbits exist
            scenes_exist.append(scene)
            if not skip_exist:
                [os.remove(os.path.join(basedir, filename)) for filename in filenames]
        #print ('scenes_exist', scenes_exist)
        
        # filter existing scenes
        if skip_exist:
            scenes = [scene for scene in scenes if scene not in scenes_exist]
            if len(scenes) == 0:
                return

        def get_urls(scenes):
            outs = []
            for scene in (scenes if isinstance(scenes, list) else [scenes]):
                out = self.template_url.format(satellite=scene[2:3], scene=scene)
                outs.append(out)
            return outs if isinstance(scenes, list) else outs[0]

        def get_patterns(subswaths, polarization, mission):
            outs = []
            for subswath in str(subswaths):
                pattern = self.template_safe.format(mission=mission.lower(),
                                                    subswath=subswath,
                                                    polarization=polarization.lower())
                outs.append(pattern)
            return outs

        def extract(url, patterns, basedir, session):
            with asf_search.remotezip(url, session) as remotezip:
                filenames = remotezip.namelist()
                for pattern in patterns:
                    #print ('pattern', pattern)
                    matching = [filename for filename in filenames if fnmatch.fnmatch(filename, pattern)]
                    for filename in matching:
                        filesize = remotezip.getinfo(filename).file_size
                        #print (filename, filesize)
                        fullname = os.path.join(basedir, filename)
                        if os.path.exists(fullname) and os.path.getsize(fullname) == filesize:
                            pass
                        else:
                            #with open(fullname, 'wb') as file:
                            #    file.write(remotezip.read(filename))
                            remotezip.extract(filename, basedir)

        if not os.path.exists(basedir):
            os.makedirs(basedir)

        urls = get_urls(scenes)
        #print ('urls', urls)

        patterns = get_patterns(subswaths, polarization, mission)
        #print ('patterns', patterns)

        # prepare authorized connection
        session = asf_search.ASFSession().auth_with_creds(self.username, self.password)

        # download scenes
        with self.tqdm_joblib(tqdm(desc='ASF Downloading Sentinel-1 SLC', total=len(urls))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(extract)(url, patterns, basedir, session) for url in urls)

        # parse scenes and convert to dataframe
        scenes = pd.DataFrame(scenes, columns=['scene'])
        scenes['time_start'] = scenes['scene'].apply(lambda name: datetime.strptime(name.split('_')[5], '%Y%m%dT%H%M%S'))
        scenes['time_end']   = scenes['scene'].apply(lambda name: datetime.strptime(name.split('_')[6], '%Y%m%dT%H%M%S'))
        scenes['mission']    = scenes['scene'].apply(lambda name: name[:3])
        scenes['orbit']      = None

        # download orbits
        if len(scenes[scenes['orbit'].isna()]):
            scenes = self._get_orbits(scenes, self.poeorb_url, session)
        if len(scenes[scenes['orbit'].isna()]):
            scenes = self._get_orbits(scenes, self.resorb_url, session)
        if len(scenes[scenes['orbit'].isna()]):
            print ('ERROR: some orbits not found')
        # download missed orbits
        asf_search.download_urls(urls=scenes.orbit.tolist(), path=basedir, session=session, processes=n_jobs)

        # return the results in a user-friendly dataframe
        scenes['scene'] = scenes.scene.apply(lambda name: self.template_url.format(satellite=name[2:3], scene=name))
        return scenes[['scene', 'orbit']]
