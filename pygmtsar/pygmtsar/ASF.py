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

    template_url = 'https://datapool.asf.alaska.edu/SLC/S{satellite}/{scene}.zip'
    template_safe = '*.SAFE/*/s1{mission}-iw{subswath}-slc-{polarization}-*'
    session = None
    
    def __init__(self, username=None, password=None):
        import asf_search
        import getpass
        if username is None:
            username = getpass.getpass('Please enter your ASF username and press Enter key:')
        if password is None:
            password = getpass.getpass('Please enter your ASF password and press Enter key:')
        self.session = asf_search.ASFSession().auth_with_creds(username, password)

    def download(self, basedir, scenes, subswaths, polarization='VV', mission='?', n_jobs=-1):
        import asf_search
        import fnmatch
        import joblib
        from tqdm.auto import tqdm
        import os

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

        def extract(url, patterns, basedir):
            with asf_search.remotezip(url, self.session) as remotezip:
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

        with self.tqdm_joblib(tqdm(desc='ASF Downloading Sentinel-1 SLC', total=len(urls))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(extract)(url, patterns, basedir) for url in urls)
