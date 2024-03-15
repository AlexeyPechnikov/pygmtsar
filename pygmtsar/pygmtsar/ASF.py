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
            # check scenes folders 
            scenes_missed = np.unique([scene for scene in scenes if f'{scene}.SAFE' not in files])
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
                            #with open(fullname, 'wb') as file:
                            #    file.write(remotezip.read(filename))
                            remotezip.extract(filename, basedir)

        # prepare authorized connection
        if session is None:
            session = self._get_asf_session()

        # download scenes
        with self.tqdm_joblib(tqdm(desc='ASF Downloading Sentinel-1 SLC', total=len(scenes_missed))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(download_scene)\
                                    (scene, subswaths, polarization, basedir, session) for scene in scenes_missed)

        # parse processed scenes and convert to dataframe
        #print ('scenes', len(scenes))
        scenes_downloaded = pd.DataFrame(scenes_missed, columns=['scene'])
        # return the results in a user-friendly dataframe
        #scenes_downloaded['scene'] = scenes_downloaded.scene\
        #                             .apply(lambda name: self.template_url.format(satellite=name[2:3], scene=name))
        return scenes_downloaded

    """
    find . -type f -name '*.tiff' -exec basename {} .tiff \; \
        | awk -F'-' -v OFS='_' '{print toupper($1), "IW_SLC__1SDV", toupper($5), toupper($6), $7, toupper($8), $9"X.SAFE"}' \
        | xargs -I {} mkdir -p {}
    """
    def download_orbits(self, basedir: str, session=None, n_jobs: int = 8, skip_exist: bool = True):
        import pandas as pd
        import requests
        import os
        import re
        import glob
        from datetime import datetime
        import asf_search
        import joblib
        from tqdm.auto import tqdm

        # create the directory if needed
        os.makedirs(basedir, exist_ok=True)

        if not skip_exist:
            orbits = glob.glob('*.EOF', root_dir=basedir)
            #print ('orbits', orbits)
            for orbit in orbits:
                os.remove(os.path.join(basedir, orbit))

        scenes = S1.scan_slc(basedir)
        scenes = scenes[scenes.orbitpath.isnull()]\
            .reset_index()\
            .groupby(['date', 'mission'])\
            .first()\
            .reset_index()[['mission', 'datetime']]
        #print ('scenes', scenes)
        if len(scenes) == 0:
            return

        # prepare authorized connection
        if session is None:
            session = self._get_asf_session()

        # download orbits
        orbits_found = []
        for product_type in ['POEORB', 'RESORB']:
            # nothing to do
            if scenes is None or len(scenes) == 0:
                continue

            url = self.resorb_url if product_type == 'RESORB' else self.poeorb_url

            # get orbits HTML list
            with tqdm(desc=f'Downloading {product_type} catalog:', total=1) as pbar:
                response = session.get(url)
                pbar.update(1)
            response.raise_for_status()
            lines = response.text.splitlines()
            orbits = [re.search(self.pattern_orbit, line).group(0) if re.search(self.pattern_orbit, line) else None for line in lines]
            orbits = pd.DataFrame(orbits, columns=['orbit']).dropna()
            orbits['mission']    = orbits['orbit'].apply(lambda name: name[:3])
            orbits['time']       = orbits['orbit'].apply(lambda name: datetime.strptime(name.split('_')[5], '%Y%m%dT%H%M%S'))
            if product_type == 'RESORB':
                orbits['time_start'] = orbits['orbit'].apply(lambda name: datetime.strptime(name.split('_')[6][1:],  '%Y%m%dT%H%M%S')) + self.offset_start
                orbits['time_end']   = orbits['orbit'].apply(lambda name: datetime.strptime(name[:-4].split('_')[7], '%Y%m%dT%H%M%S')) - self.offset_end
            elif product_type == 'POEORB':
                # use daily orbits
                orbits['time_start'] = orbits['orbit'].apply(lambda name: datetime.strptime(name.split('_')[6][1:][:9]  + '235959', '%Y%m%dT%H%M%S'))
                orbits['time_end']   = orbits['orbit'].apply(lambda name: datetime.strptime(name[:-4].split('_')[7][:9] + '000001', '%Y%m%dT%H%M%S'))
            #print ('orbits', orbits.head(3))

            urls = []
            scenes_missed = []
            for scene in scenes.itertuples():
                #print ('scene', scene)
                # look for the orbit file
                orbit = orbits[(orbits.mission == scene.mission)&\
                               (orbits.time_start <= scene.datetime)&\
                               (orbits.time_end   >= scene.datetime)].sort_values('time')
                #print ('orbit', orbit)
                # add the recent orbit
                orbit = orbit.tail(1).orbit.item() if len(orbit) >= 1 else None
                if orbit is None:
                    scenes_missed.append(scene)
                else:
                    urls.append(url + orbit)
                    orbits_found.append(orbit)

            # this routine can use multiple threads but it does not provide a progress indicator
            with self.tqdm_joblib(tqdm(desc=f'ASF Downloading {product_type} Orbits', total=len(orbits_found))) as progress_bar:
                joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(asf_search.download_urls)(urls=[url], path=basedir, session=session) \
                                               for url in urls)

            # remove processed scenes
            if len(scenes_missed) > 0:
                scenes = pd.DataFrame(scenes_missed)[['mission', 'datetime']]
            else:
                scenes = None

        # check for scenes without orbits
        if scenes is not None and len(scenes) > 0:
            raise Exception(f'ERROR: missed {len(scenes)} orbits')
        return pd.Series(orbits_found)

    def download(self, *args, **kwarg):
        print ('NOTE: Function is deprecated. Use ASF.download_scenes() and ASF.download_orbits().')
