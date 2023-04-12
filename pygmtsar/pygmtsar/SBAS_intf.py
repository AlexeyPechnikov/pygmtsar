#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_topo_ra import SBAS_topo_ra

class SBAS_intf(SBAS_topo_ra):

    def intf(self, subswath, pair, **kwargs):
        import os

        # extract dates from pair
        date1, date2 = pair

        prm_ref = self.PRM(subswath, date1)
        prm_rep = self.PRM(subswath, date2)

        topo_ra_file = os.path.join(self.basedir, f'F{subswath}_topo_ra.grd')
        #print ('SBAS intf kwargs', kwargs)
        prm_ref.intf(prm_rep,
                     basedir=self.basedir,
                     topo_ra_fromfile = topo_ra_file,
                     **kwargs)

    def intf_parallel(self, pairs, n_jobs=-1, **kwargs):
        import pandas as pd
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        import os

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        subswaths = self.get_subswaths()

        # for now (Python 3.10.10 on MacOS) joblib loads the code from disk instead of copying it
        kwargs['chunksize'] = self.chunksize
        
        # this way does not work properly for long interferogram series on MacOS
        # see https://github.com/mobigroup/gmtsar/commit/3eea6a52ddc608639e5e06306bce2f973a184fd6
        #with self.tqdm_joblib(tqdm(desc='Interferograms', total=len(pairs)*len(subswaths))) as progress_bar:
        #    joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.intf)(subswath, pair, **kwargs) \
        #        for subswath in subswaths for pair in pairs)

        # workaround: start a set of jobs together but not more than available cpu cores at once
        from joblib.externals import loky
        if n_jobs == -1:
            n_jobs = joblib.cpu_count()
        # create list of arrays [subswath, date1, date2] where all the items are strings
        subpairs = [[subswath, pair[0], pair[1]] for subswath in subswaths for pair in pairs]
        n_chunks = int(np.ceil(len(subpairs)/n_jobs))
        chunks = np.array_split(subpairs, n_chunks)
        #print ('n_jobs', n_jobs, 'n_chunks', n_chunks, 'chunks', [len(chunk) for chunk in chunks])
        with tqdm(desc='Interferograms', total=len(subpairs)) as pbar:
            for chunk in chunks:
                loky.get_reusable_executor(kill_workers=True).shutdown(wait=True)
                with joblib.parallel_backend('loky', n_jobs=n_jobs):
                    # convert string subswath to integer value
                    joblib.Parallel()(joblib.delayed(self.intf)(int(subswath), [date1, date2], **kwargs) \
                        for (subswath,date1,date2) in chunk)
                    pbar.update(len(chunk))
