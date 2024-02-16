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
    template_safe = '*.SAFE/*/{mission}-iw{subswath}-slc-{polarization}-*'
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
    resorb_start_offset = timedelta(seconds=(12 * 86400.0) // 175.0 + 60)
    resorb_end_offset = timedelta(seconds=300)
    poeorb_start_offset = timedelta(seconds=(12 * 86400.0) // 175.0 + 60)
    poeorb_end_offset = timedelta(seconds=300)
    
    def __init__(self, username=None, password=None):
        import asf_search
        import getpass
        if username is None:
            username = getpass.getpass('Please enter your ASF username and press Enter key:')
        if password is None:
            password = getpass.getpass('Please enter your ASF password and press Enter key:')
        self.username = username
        self.password = password

    def _get_orbits(self, scenes, session, kind):
        import pandas as pd
        import re
        from datetime import datetime

        if kind == 'RESORB':
            url = self.resorb_url
            orbit_start_offset = self.resorb_start_offset
            orbit_end_offset = self.resorb_end_offset
        elif kind == 'POEORB':
            url = self.poeorb_url
            orbit_start_offset = self.poeorb_start_offset
            orbit_end_offset = self.poeorb_end_offset
        else:
            raise Exception(f'Invalid orbit type: {kind}.')
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
            #print ('scene', scene)
            if scene.orbit is not None:
                continue
            # look for the orbit file
            orbit = orbits[(orbits.mission == scene.mission)&\
                           (orbits.time_start <= scene.time_start - orbit_start_offset)&\
                           (orbits.time_end   >= scene.time_end   + orbit_end_offset  )].sort_values('time')
            #print ('orbit', orbit)
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
    def download(self, basedir, scenes, subswaths, polarization='VV', n_jobs=4, skip_exist=True):
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

        # collect all the existing files once
        files = os.listdir(basedir)
        #print ('files', len(files))

        # remove all orbits and re-download the same or new ones
        if not skip_exist:
            orbits = [res.group(0) if res else None for file in files for res in [re.search(self.pattern_orbit, file)]]
            orbits = list(filter(None, orbits))
            for orbit in orbits:
                os.remove(os.path.join(basedir, orbit))

        # detect all the local orbits for all scenes
        scenes_orbits = []
        # skip existing scenes
        if skip_exist:
            # output dataframe including scenes and orbits so the both required
            # missed scene or orbit means scene as missed
            scenes_missed = []
            for scene in scenes:
                # check scene
                if not f'{scene}.SAFE' in files:
                    scenes_missed.append(scene)

                # check orbit
                mission = scene[:3]
                dt = datetime.strptime(scene.split('_')[5], '%Y%m%dT%H%M%S')
                day        = dt.strftime('%Y%m%d')
                day_before = (dt + timedelta(days=-1)).strftime('%Y%m%d')
                day_after  = (dt + timedelta(days=1)).strftime('%Y%m%d')
                # check presision orbits
                pattern    = self.template_orbit.format(mission=mission, date_start=day_before, date_end=day_after)
                #print ('Precision orbit pattern', pattern)
                filenames  = sorted(filter(re.compile(pattern).match, files))
                [scenes_orbits.append((scene, filename)) for filename in filenames]
                #print ('filenames', filenames)
                # check approximate orbits when start time is right after midnight and the start date is day before
                if len(filenames) == 0:
                    #print ('Precision orbit not found')
                    pattern   = self.template_orbit.format(mission=mission, date_start=day_before, date_end=day)
                    #print ('Approximate orbit pattern', pattern)
                    filenames = sorted(filter(re.compile(pattern).match, files))
                    [scenes_orbits.append((scene, filename)) for filename in filenames]
                # check approximate orbits
                if len(filenames) == 0:
                    #print ('Precision orbit not found')
                    pattern   = self.template_orbit.format(mission=mission, date_start=day, date_end=day)
                    #print ('Approximate orbit pattern', pattern)
                    filenames = sorted(filter(re.compile(pattern).match, files))
                    [scenes_orbits.append((scene, filename)) for filename in filenames]
                # orbit file is missed when presision and approximate orbits not found
                if len(filenames) == 0:
                    #print ('Approximate orbit not found')
                    scenes_missed.append(scene)

            #print ('scenes_missed', len(scenes_missed))
            scenes = np.unique(scenes_missed)
            #print ('scenes_orbits', scenes_orbits)

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
                            #with open(fullname, 'wb') as file:
                            #    file.write(remotezip.read(filename))
                            remotezip.extract(filename, basedir)

        # do not initialize internet connection, it allows to work offline when all the scenes and orbits already available
        if len(scenes) == 0:
            return

        # prepare authorized connection
        session = asf_search.ASFSession().auth_with_creds(self.username, self.password)

        # download scenes
        with self.tqdm_joblib(tqdm(desc='ASF Downloading Sentinel-1 SLC', total=len(scenes))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(download_scene)(scene, subswaths, polarization, basedir, session) for scene in scenes)
        # debug
        #[download_scene(scene, subswaths, polarization, basedir, session) for scene in scenes]

        # parse processed scenes and convert to dataframe
        #print ('scenes', len(scenes))
        scenes = pd.DataFrame(scenes, columns=['scene'])
        scenes['time_start'] = scenes['scene'].apply(lambda name: datetime.strptime(name.split('_')[5], '%Y%m%dT%H%M%S'))
        scenes['time_end']   = scenes['scene'].apply(lambda name: datetime.strptime(name.split('_')[6], '%Y%m%dT%H%M%S'))
        scenes['mission']    = scenes['scene'].apply(lambda name: name[:3])
        scenes['orbit']      = None

        # download orbit catalogs and look for the required scenes orbits
        if len(scenes[scenes['orbit'].isna()]):
            with tqdm(desc='Downloading POEORB catalog:', total=1) as pbar:
                scenes = self._get_orbits(scenes, session, 'POEORB')
                pbar.update(1)
        if len(scenes[scenes['orbit'].isna()]):
            with tqdm(desc='Downloading RESORB catalog:', total=1) as pbar:
                scenes = self._get_orbits(scenes, session, 'RESORB')
                pbar.update(1)
        if len(scenes[scenes['orbit'].isna()]):
            print ('ERROR: some orbits not found locally and in orbit catalogs')

        def download_orbit(scene, url, basedir, session):
            #print ('download_orbit', scene, url)
            # cleanup one or multiple existing orbits
            # orbits can be outdated and have different names
            # detect all the existing orbits for the scene
            orbits = [scene_orbit[1] for scene_orbit in scenes_orbits if scene_orbit[0]==scene]       
            fullnames = [os.path.join(basedir, orbit) for orbit in orbits]
            #print ('fullnames', fullnames)
            # it can be that scene is missed but the orbit esists when skip_exist=True
            for fullname in fullnames:
                if os.path.exists(fullname):
                    os.remove(fullname)
            # download the same or a new orbit 
            # this routine can use multiple threads but it does not provide a progress indicator
            asf_search.download_urls(urls=[url], path=basedir, session=session)

        # download orbits
        with self.tqdm_joblib(tqdm(desc='ASF Downloading Sentinel-1 Orbits', total=len(scenes.orbit))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(download_orbit)(scene, url, basedir, session) \
                                           for (scene, url) in zip(scenes.scene, scenes.orbit))
        # debug
    #    [download_orbit(url, basedir, session) for url in scenes.orbit]

        # return the results in a user-friendly dataframe
        scenes['scene'] = scenes.scene.apply(lambda name: self.template_url.format(satellite=name[2:3], scene=name))
        return scenes[['scene', 'orbit']]
