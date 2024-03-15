# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2024, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .tqdm_joblib import tqdm_joblib
from .S1 import S1

class ESA(tqdm_joblib):
    from datetime import timedelta, datetime

    # Default URL endpoint for performing user authentication with CDSE
    token_url = 'https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token'
    template_query = (
        "startswith(Name,'{mission_id}') and contains(Name,'{orbit_type}') "
        "and ContentDate/Start lt '{start_time}' and ContentDate/End gt '{stop_time}'"
    )
    # Default URL endpoint for the Copernicus Data Space Ecosystem (CDSE) query REST service
    query_url = 'https://catalogue.dataspace.copernicus.eu/odata/v1/Products'
    # Default URL endpoint for CDSE download REST service
    download_url = 'https://zipper.dataspace.copernicus.eu/odata/v1/Products'
    #pattern_orbit  = r'S1\w_OPER_AUX_\w{3}ORB_OPOD_\d{8}T\d{6}_V\d{8}T\d{6}_\d{8}T\d{6}.EOF'
    #pattern_orbit  = r'S1\w_OPER_AUX_(POE|RES)ORB_OPOD_\d{8}T\d{6}_V\d{8}T\d{6}_\d{8}T\d{6}.EOF'
    #template_orbit = '{mission}_OPER_AUX_(POE|RES)ORB_OPOD_\d+T\d+_V{date_start}T\d+_{date_end}T\d+.EOF'
    # see _select_orbit.py in sentineleof package
    #Orbital period of Sentinel-1 in seconds
    #T_ORBIT = (12 * 86400.0) / 175.0
    offset_start = timedelta(seconds=(12 * 86400.0) // 175.0 + 60)
    offset_end = timedelta(seconds=300)
    
    def __init__(self, username: str = None, password: str= None):
        """
        Register your ESA Copernicus datastore account at https://dataspace.copernicus.eu/
        """
        import asf_search
        import getpass
        if username is None:
            username = getpass.getpass('Please enter your Copernicus Data Space username and press Enter key:')
        if password is None:
            password = getpass.getpass('Please enter your Copernicus Data Space password and press Enter key:')
        self.username = username
        self.password = password

    def _get_copernicus_access_token(self) -> str:
        """
        Authenticate with the Copernicus Data Space Identity Server to obtain an access token.

        Parameters:
        - username: Copernicus Data Space username.
        - password: Copernicus Data Space password.

        Returns:
        - A string containing the access token.
        
        Example:
        username = 'your_copernicus_username'
        password = 'your_copernicus_password'
        access_token = get_copernicus_access_token(username, password)
        print("Access Token:", access_token)
        """
        import requests
        
        # Endpoint for obtaining the access token
        data = {
            'client_id': 'cdse-public',
            'username': self.username,
            'password': self.password,
            'grant_type': 'password',
        }
        try:
            response = requests.post(self.token_url, data=data)
            response.raise_for_status()
            access_token = response.json().get('access_token')
            if not access_token:
                raise Exception('Failed to obtain access token.')
            return access_token
        except requests.exceptions.RequestException as e:
            raise Exception(f'Authentication failed: {e}')

    def _download_orbit(self, access_token : str, time: datetime,
                       mission_id: str, basedir: str) -> bool:
        """
        Queries and downloads RESORB/POEORB files for the given satellite within the specified time range.
        """
        import requests
        import os
        
        for product_type in ['AUX_POEORB', 'AUX_RESORB']:
            query = self.template_query.format(
                start_time = (time - self.offset_start).strftime('%Y-%m-%dT%H:%M:%SZ') if product_type == 'AUX_RESORB' else
                         time.replace(hour=0,  minute=0,  second=1 ).strftime('%Y-%m-%dT%H:%M:%SZ'),
                stop_time  = (time + self.offset_end).strftime('%Y-%m-%dT%H:%M:%SZ')   if product_type == 'AUX_RESORB' else
                         time.replace(hour=23, minute=59, second=59).strftime('%Y-%m-%dT%H:%M:%SZ'),
                mission_id = mission_id,
                orbit_type = product_type,
            )
            #print ('query', query)
        
            # Query the Copernicus Data Space for available orbit files
            headers = {'Authorization': f'Bearer {access_token}'}
            response = requests.get(self.query_url, params={'$filter': query}, headers=headers)
            response.raise_for_status()

            # Process the query results
            results = response.json().get('value', [])
            if not results:
                print(f'No {product_type} files found to cover time {time}')
                continue
        
            if len(results) > 1:
                print(f'Found multiple {product_type} orbits for time {time}', results)
            # TODO: what is the sorting for multiple orbits?
            result = results[0]
            
            # Download the orbit file
            orbit_file_id = result['Id']
            orbit_file_name = result['Name']
            download_url = f'{self.download_url}({orbit_file_id})/$value'
            response = requests.get(download_url, headers=headers, stream=True)
            response.raise_for_status()
            file_path = os.path.join(basedir, orbit_file_name)
#             with open(file_path, 'wb') as f:
#                 for chunk in response.iter_content(chunk_size=256*1024):
#                     f.write(chunk)
            with open(file_path, 'wb') as f:
                f.write(response.content)
            return result['Name']

        return

    def download_orbits(self, basedir: str, n_jobs: int = 2, skip_exist: bool = True, debug: bool = False):
        import pandas as pd
        import os
        import glob
        from tqdm.auto import tqdm
        import joblib

        # create the directory if needed
        os.makedirs(basedir, exist_ok=True)

        if not skip_exist:
            orbits = glob.glob('*.EOF', root_dir=basedir)
            #print ('orbits', orbits)
            for orbit in orbits:
                os.remove(os.path.join(basedir, orbit))

        scenes = S1.scan_slc(basedir)
        mission_date = scenes[scenes.orbitpath.isnull()]\
            .reset_index()\
            .groupby(['date', 'mission'])\
            .first()\
            .reset_index()[['mission', 'datetime']]
        #print ('mission_date', mission_date)
        if len(mission_date) == 0:
            return

        access_token = self._get_copernicus_access_token()
        # download orbits
        joblib_backend = 'sequential' if debug else 'loky'
        with self.tqdm_joblib(tqdm(desc='ESA Downloading Sentinel-1 Orbits', total=len(mission_date))) as progress_bar:
            orbits = joblib.Parallel(n_jobs=n_jobs, backend=joblib_backend)(joblib.delayed(self._download_orbit)\
                                    (access_token, date, mission, basedir)
                                           for (mission, date) in mission_date.values)
        return pd.Series(orbits)
