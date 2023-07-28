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

    def get_adi(self, subswath=None, threshold=None, chunksize=None):
        """
        TODO: threshold
        """
    
        adis = self.open_grids(None, 'adi', chunksize=chunksize)
        return adis[0] if len(adis)==1 else adis
    
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
    def adi_parallel(self, dates=None, chunksize=None):
        import xarray as xr
        import numpy as np
        import dask
        import os

        if chunksize is None:
            chunksize = self.chunksize

        for subswath in self.get_subswaths():
            if dates is None:
                dates = self.df[self.df['subswath']==subswath].index.values
            #print ('dates', dates)
            # select radar coordinates extent
            rng_max = self.PRM(subswath).get('num_rng_bins')
            #print ('azi_max', azi_max, 'rng_max', rng_max)
            # use SLC-related chunks for faster processing
            minichunksize = int(np.round(chunksize**2/rng_max))
            #print (minichunksize)
            amps = [self.PRM(subswath, date).read_SLC_int(amplitude=True, chunksize=minichunksize) for date in dates]
            # build stack
            amps = xr.concat(amps, dim='date')
            # normalize image amplitudes
            #mean = amps.mean(dim=['y','x']).compute()
            tqdm_dask(mean := dask.persist(amps.mean(dim=['y','x'])), desc=f'Amplitude Normalization sw{subswath}')
            #print ('mean', mean)
            # dask.persist returns tuple
            norm = (mean[0]/mean[0].mean(dim='date'))
            del mean
            #print ('norm', norm)
            # compute Amplitude Dispersion Index (ADI)
            #adi = ((norm*amps).std(dim='date')/(norm*amps).mean(dim='date')).rename('adi')
            #del norm, amps
            stats = (norm*amps).pipe(lambda x: (x.mean(dim='date'), x.std(dim='date')))
            del amps, norm
            adi = (stats[1] / stats[0]).rename('adi')
            del stats
            # Define the encoding and choose the appropriate storage scheme
            adi_filename = self.get_filenames(subswath, None, f'adi')
            if os.path.exists(adi_filename):
                os.remove(adi_filename)
            # save by fastest way using SLC-related chunks like (48, 1024)
            delayed = adi.rename({'y': 'a', 'x': 'r'}).to_netcdf(adi_filename,
                         engine=self.engine,
                         encoding={'adi': self.compression(adi.shape, chunksize=(minichunksize, chunksize))},
                         compute=False)
            tqdm_dask(dask.persist(delayed), desc=f'Amplitude Dispersion Index (ADI) sw{subswath}')
            del adi, delayed
            # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
            import gc; gc.collect()
