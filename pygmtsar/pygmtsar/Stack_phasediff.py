# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_topo import Stack_topo
from .tqdm_dask import tqdm_dask

class Stack_phasediff(Stack_topo):

    def phasediff(self, pair, topo_fromfile=None, chunksize=None, debug=False):
        import os
        import pandas as pd
        import numpy as np
        import xarray as xr
        import dask.array
        import warnings
        from datetime import datetime
        # suppress Dask warning "RuntimeWarning: divide by zero encountered in divide"
        #warnings.filterwarnings("ignore", category=RuntimeWarning, module="dask.core")

        # unique filenames specifier
        #timenow = datetime.now().strftime("%F_%T.%f").replace(':', '.')

        # define lost class variables due to joblib
        if chunksize is None:
            chunksize = self.chunksize

        # convert to 2D single-element array
        pairs = self.get_pairs([pair] if np.asarray(pair).ndim == 1 else pair)
        pair = pairs[['ref','rep']].astype(str).values[0]

        # extract dates from pair
        date1, date2 = pair

        prm_ref = self.PRM(date1)
        prm_rep = self.PRM(date2)
        if debug:
            print ('PRM reference:', prm_ref.filename)
            print ('PRM repeat:   ', prm_rep.filename)

        if topo_fromfile is None:
            topo_fromfile = self.get_filename('topo')

        subswath = self.get_subswath()

        # make full file name, use workaround for 'weight' argument name defined without extension
        fullname = lambda name: os.path.join(self.basedir, f'F{subswath}_{date1}_{date2}_{name}'.replace('-',''))

        # prepare PRMs for the calculation below
        if debug:
            print ('SAT_baseline:\n', prm_ref.SAT_baseline(prm_rep, tail=9))
        prm_rep.set(prm_ref.SAT_baseline(prm_rep, tail=9))
        prm_ref.set(prm_ref.SAT_baseline(prm_ref).sel('SC_height','SC_height_start','SC_height_end'))

        # for topo use relative path from PRM files directory
        # use imag.grd=bf for GMT native, C-binary format
        # cleanup
        for name in ['real.grd', 'imag.grd']:
            filename = fullname(name)
            if os.path.exists(filename):
                os.remove(filename)
        prm_ref.phasediff(prm_rep, topo_fromfile=topo_fromfile,
                       imag_tofile=fullname('imag.grd'),
                       real_tofile=fullname('real.grd'),
                       debug=debug)

    #     # original SLC (do not flip vertically)
    #     amp1 = prm_ref.read_SLC_int(intensity=True)
    #     amp2 = prm_rep.read_SLC_int(intensity=True)
    #     # phasediff tool output files (flip vertically)
    #     imag = xr.open_dataarray(fullname('imag.grd'), engine=self.engine, chunks=chunksize)
    #     imag.data = dask.array.flipud(imag)
    #     real = xr.open_dataarray(fullname('real.grd'), engine=self.engine, chunks=chunksize)
    #     real.data = dask.array.flipud(real)

    #     out = xr.Dataset({'amp1': amp1, 'amp2': amp2, 'phasediff': real +1j*imag})
    #     out['ref'] = date1
    #     out['rep'] = date2
    #     out['pair'] = f'{date1} {date2}'
    #     del amp1, amp2, real, imag
    #     return out

    def stack_phasediff(self, pairs, topo_fromfile=None, chunksize=None, n_jobs=-1):
        """
        Build phase difference for all pairs.

        Parameters
        ----------
        pairs : list
            List of date pairs (baseline pairs).
        topo_fromfile : str, optional
            Define radar coordinate topography. Default is topography created by geocode().
        n_jobs : int, optional
            Number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
        Returns
        -------
        None
        """
        from tqdm.auto import tqdm
        import joblib

        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs = self.get_pairs(pairs)[['ref', 'rep']].astype(str).values

        with self.tqdm_joblib(tqdm(desc='Phase Difference', total=len(pairs))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.phasediff)(pair, topo_fromfile, chunksize) \
                for pair in pairs)
