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
        from joblib.externals import loky
        import os

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        subswaths = self.get_subswaths()

        # this way does not work properly for long interferogram series
        #with self.tqdm_joblib(tqdm(desc='Interferograms', total=len(pairs))) as progress_bar:
        #    joblib.Parallel(n_jobs=-1)(joblib.delayed(self.intf)(pair, **kwargs) for pair in pairs)

        # start a set of jobs together but not more than available cpu cores at once
        if n_jobs == -1:
            n_jobs = joblib.cpu_count()
        n_chunks = int(np.ceil(len(pairs)/n_jobs))
        chunks = np.array_split(pairs, n_chunks)
        #print ('n_jobs', n_jobs, 'n_chunks', n_chunks, 'chunks', [len(chunk) for chunk in chunks])
        with tqdm(desc='Interferograms', total=len(pairs)*len(subswaths)) as pbar:
            for chunk in chunks:
                loky.get_reusable_executor(kill_workers=True).shutdown(wait=True)
                with joblib.parallel_backend('loky', n_jobs=n_jobs, inner_max_num_threads=1):
                    joblib.Parallel()(joblib.delayed(self.intf)(subswath, pair, **kwargs) \
                        for subswath in subswaths for pair in chunk)
                pbar.update(len(chunk)*len(subswaths))

        # backward compatibility wrapper
        # for a single subswath don't need to call SBAS.merge_parallel()
        # for subswaths merging and total coordinate transformation matrices creation 
        #if len(subswaths) == 1:
        #    # build geo transform matrices for interferograms
        #    self.transforms(subswaths[0], pairs)
