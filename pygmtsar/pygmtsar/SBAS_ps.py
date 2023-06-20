# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_stl import SBAS_stl
from .tqdm_dask import tqdm_dask

class SBAS_ps(SBAS_stl):
    
    # define PS candidates using Amplitude Dispersion Index (ADI)
    def ps_parallel(self, dates=None, threshold=0.25, chunksize=None, debug=False):
        import xarray as xr
        import numpy as np

        if dates is None:
            dates = self.df.index

        # scale factor is selected to have the same number of pixels per square chunk
        if chunksize is None:
            chunksize = int(np.round(self.chunksize**2/self.PRM().get('num_rng_bins')))

        amps = []
        for subswath in self.get_subswaths():
            amps_subswath = []
            for date in dates:
                #print (subswath, date)
                prm = self.PRM(subswath, date)
                amp = prm.read_SLC_int(amplitude=True, chunksize=chunksize)
                amp['date'] = date
                amps_subswath.append(amp)
            # build stack
            amps_subswath = xr.concat(amps_subswath, dim='date')
            # compute Amplitude Dispersion Index (ADI)
            adi = amps_subswath.std(dim='date')/amps_subswath.mean(dim='date')
            # use the specified threshold
            amps.append(adi.where(adi<=threshold))

        return amps[0] if len(amps)==1 else amps
