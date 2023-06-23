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

    #from pygmtsar import tqdm_dask
    #SBAS.ps_parallel = ps_parallel    
    #sbas.ps_parallel(interactive=True)
    #sbas.ps_parallel()
    #adi = sbas.open_grids(None, 'ps')
    #adi
    #ps_decimator = sbas.pixel_decimator(resolution_meters=60, grid=adi, debug=True)
    #adi_dec = adi.coarsen({'y': 4, 'x': 16}, boundary='trim').min()
    #adi_dec
    # define PS candidates using Amplitude Dispersion Index (ADI)
    def ps_parallel(self, dates=None, threshold=None, chunksize=None, interactive=False, debug=False):
        import xarray as xr
        import numpy as np
        import dask
        import os

        if dates is None:
            dates = self.df.index

        # scale factor is selected to have the same number of pixels per square chunk
        if chunksize is None:
            chunksize = self.chunksize
        # use SLC-related chunks for faster processing
        bigchunk = self.PRM().get('num_rng_bins')
        minichunksize = int(np.round(chunksize**2/bigchunk))
        #print (minichunksize, bigchunk)

        amps = []
        for subswath in self.get_subswaths():
            subamps = []
            for date in dates:
                #print (subswath, date)
                prm = self.PRM(subswath, date)
                amp = prm.read_SLC_int(amplitude=True, chunksize=minichunksize)
                amp['date'] = date
                subamps.append(amp)
            # build stack
            subamps = xr.concat(subamps, dim='date')
            # normalize image amplitudes
            mean = subamps.mean(dim=['y','x'])
            norm = mean/mean.mean(dim='date')
            # compute Amplitude Dispersion Index (ADI)
            adi = (norm*subamps).std(dim='date')/(norm*subamps).mean(dim='date')
            # use the specified threshold
            if threshold is not None:
                amps.append(adi.where(adi<=threshold))
            else:
                amps.append(adi)

        if interactive:
            return amps[0] if len(amps)==1 else amps

        delayeds = []
        for (subswath, ps) in zip(self.get_subswaths(), amps):
            # Define the encoding and choose the appropriate storage scheme
            ps_filename = self.get_filenames(subswath, None, f'ps')
            if os.path.exists(ps_filename):
                os.remove(ps_filename)
            # save by fastest way using chunks like to (12, 512)
            delayed = ps.rename('adi').to_netcdf(ps_filename,
                         engine=self.engine,
                         encoding={'adi': self.compression(ps.shape, chunksize=(minichunksize, 512))},
                         compute=False)
            delayeds.append(delayed)
        tqdm_dask(dask.persist(delayeds), desc='Persistent Scatterers')
