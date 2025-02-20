# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .tqdm_joblib import tqdm_joblib

class EOF(tqdm_joblib):

    import pandas as pd
    from datetime import timedelta

    http_timeout = 30
    #https://s1orbits.pechnikov.workers.dev/S1A/2014/04/07/index.csv
    #https://s1orbits.pechnikov.workers.dev/S1A/2014/04/07/S1A_OPER_AUX_POEORB_OPOD_20210301T130653_V20140406T225944_20140408T005944.EOF.zip
    orbits_url = 'https://s1orbits.pechnikov.workers.dev/{mission}/{year}/{month:02}/{day:02}/'
    # see _select_orbit.py in sentineleof package
    #Orbital period of Sentinel-1 in seconds
    #T_ORBIT = (12 * 86400.0) / 175.0
    # ESA orbits doe not follow the specification and some orbits are missed
    # example scene: S1A_IW_SLC__1SDV_20240505T002709_20240505T002733_053728_06870F_8D20
    #orbit_offset_start = timedelta(seconds=(12 * 86400.0) // 175.0 + 60)
    # less strict rule allows to find the required orbits
    orbit_offset_start = timedelta(seconds=3600)
    orbit_offset_end = timedelta(seconds=300)

    @staticmethod
    def download(basedir: str, scenes: list | pd.DataFrame,
                        n_jobs: int = 8, joblib_backend='loky', skip_exist: bool = True,
                        retries: int = 30, timeout_second: float = 3):
        """
        Downloads orbit files corresponding to the specified Sentinel-1 scenes.
    
        Parameters
        ----------
        basedir : str
            The directory where the downloaded orbit files will be saved.
        scenes : list or pandas.DataFrame
            List of scene identifiers or a DataFrame containing scenes for which the orbits are to be downloaded.
        n_jobs : int, optional
            The number of concurrent download jobs. Default is 8.
        joblib_backend : str, optional
            The backend for parallel processing. Default is 'loky'.
        skip_exist : bool, optional
            If True, skips downloading orbits that already exist. Default is True.
    
        Returns
        -------
        pandas.Series
            A Series containing the names of the downloaded orbit files.
    
        Raises
        ------
        ValueError
            If an invalid scenes argument is provided or no suitable orbit files are found.
        """
        import pandas as pd
        import requests
        import os
        import re
        import glob
        from datetime import datetime
        import time
        import joblib
        import zipfile
        from io import BytesIO
        import xmltodict
        from tqdm.auto import tqdm
    
        # create the directory if needed
        os.makedirs(basedir, exist_ok=True)
    
        if not skip_exist:
            orbits = glob.glob('*.EOF', root_dir=basedir)
            #print ('orbits', orbits)
            for orbit in orbits:
                os.remove(os.path.join(basedir, orbit))
    
        if isinstance(scenes, pd.DataFrame):
            if skip_exist:
                # ignore scenes with orbits
                df_data = scenes[scenes.orbit.isnull()][['mission', 'datetime']].values
            else:
                # process all the scenes
                df_data = scenes[['mission', 'datetime']].values
        elif isinstance(scenes, list):
            df_data = [(scene.split('_')[0], datetime.strptime(scene.split('_')[5],'%Y%m%dT%H%M%S')) for scene in scenes]
        else:
            raise ValueError(f'Expected secenes argument is list or Pandas DataFrame')
        # nothing to do
        if len(df_data) == 0:
            return
        df = pd.DataFrame(df_data, columns=['mission', 'datetime'])
        df['date'] = df['datetime'].dt.date
        df = df.groupby(['date', 'mission'])['datetime'].first().reset_index()
        del df['date']
        
        def download_index(mission, dt, orbit_offset_start, orbit_offset_end):
            url = EOF.orbits_url.format(mission=mission,
                                       year=dt.date().year,
                                       month=dt.date().month,
                                       day=dt.date().day)
            #print ('url', url + 'index.csv')
            with requests.get(url + 'index.csv', timeout=EOF.http_timeout) as response:
                response.raise_for_status()
                # the most recent orbits preferred
                lines = sorted(response.text.splitlines(), reverse=True)
                orbits = pd.DataFrame(lines, columns=['orbit']).dropna()
                orbits['product']    = orbits['orbit'].apply(lambda name: name.split('_')[3])
                orbits['mission']    = orbits['orbit'].apply(lambda name: name[:3])
                orbits['time']       = orbits['orbit'].apply(lambda name: datetime.strptime(name.split('_')[5], '%Y%m%dT%H%M%S'))
            # detect suitable orbit and select the first one
            for orbit in orbits.itertuples():
                #print ('orbit', orbit)
                orbit_parts = orbit.orbit.split('.')[0].split('_')
                #print ('orbit_parts', orbit_parts)
                if orbit.product == 'RESORB':
                    # select the one covering the interval
                    time_start = datetime.strptime(orbit_parts[6][1:],  '%Y%m%dT%H%M%S') + orbit_offset_start
                    time_end   = datetime.strptime(orbit_parts[7], '%Y%m%dT%H%M%S') - orbit_offset_end
                    #print ('time_start', time_start, 'time_end', time_end)
                    if not ((time_start <= dt) & (time_end >= dt)):
                        continue
                    return (orbit.orbit, url)
                elif orbit.product == 'POEORB':
                    # select the most recent one
                    return (orbit.orbit, url)
                else:
                    # in case of data parse or another error
                    raise ValueError(f'Unexpected orbit product {orbit.product}')
            # downloading is not possible
            raise ValueError(f'Orbit product not found for mission {mission} and timestamp {dt}')
    
        # TODO: unzip files
        def download_orbit(basedir, url):
            filename = os.path.join(basedir, os.path.basename(os.path.splitext(url)[0]))
            #print ('url', url, 'filename', filename)
            with requests.get(url, timeout=EOF.http_timeout) as response:
                response.raise_for_status()
                with open(filename, 'wb') as f:
                    with zipfile.ZipFile(BytesIO(response.content), 'r') as zip_in:
                        zip_files = zip_in.namelist()
                        if len(zip_files) == 0:
                            raise Exception('ERROR: Downloaded file is empty zip archive.')
                        if len(zip_files) > 1:
                            raise Exception('NOTE: Downloaded zip archive includes multiple files.')
                        # extract specific file content
                        orbit_content = zip_in.read(zip_files[0])
                        # check XML validity
                        xmltodict.parse(orbit_content)
                        f.write(orbit_content)

        def download_with_retry(func, retries, timeout_second, *args, **kwargs):
            for retry in range(retries):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    print(f'ERROR: download attempt {retry+1} failed: {e}')
                    if retry + 1 >= retries:
                        raise
                time.sleep(timeout_second)

        # download orbits index files and detect the orbits.
        # joblib reads the class file from disk,
        # and to allow modification of orbit_offset_* on the fly, add them to the arguments
        with EOF.tqdm_joblib(tqdm(desc='Downloading Sentinel-1 Orbits Index:', total=len(df))) as progress_bar:
            orbits = joblib.Parallel(n_jobs=n_jobs, backend=joblib_backend)(joblib.delayed(download_with_retry)\
                                    (download_index, retries, timeout_second,
                                    mission=scene.mission,
                                    dt=scene.datetime,
                                    orbit_offset_start=EOF.orbit_offset_start,
                                    orbit_offset_end=EOF.orbit_offset_end) for scene in df.itertuples())
        # convert to dataframe for processing
        orbits = pd.DataFrame(orbits, columns=['orbit', 'url'])
        # exclude duplicates if needed 
        #orbits = orbits.groupby(['orbit']).first().reset_index()
        
        # download orbits index files and detect the orbits
        with EOF.tqdm_joblib(tqdm(desc='Downloading Sentinel-1 Orbits:', total=len(orbits))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs, backend=joblib_backend)(joblib.delayed(download_with_retry)\
                                    (download_orbit, retries, timeout_second,
                                    basedir=basedir,
                                    url=orbit.url + orbit.orbit) for orbit in orbits.itertuples())
        return orbits['orbit']
