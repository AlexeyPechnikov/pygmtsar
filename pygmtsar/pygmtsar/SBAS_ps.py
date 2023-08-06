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

    def get_adi(self, subswath=None, chunksize=None):

        if subswath is None:
            subswaths = self.get_subswaths()
        else:
            subswaths = [subswath]

        if chunksize is None:
            chunksize = self.chunksize

        adis = []
        for subswath in subswaths:
            filename = f'F{subswath}_adi'
            adi = self.open_model(filename, chunksize=chunksize)
            if 'a' in adi.dims and 'r' in adi.dims:
                adis.append(adi.rename({'a': 'y', 'r': 'x'}))
            else:
                adis.append(adi)

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
    def adi_parallel(self, dates=None, scale=1e9, chunksize=None):
        """
        scale amplitude to make minimum valid value about 1
        """
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
            amps = [self.PRM(subswath, date).read_SLC_int(intensity=True, chunksize=minichunksize) for date in dates]
            # build stack
            amps = xr.concat(amps, dim='date')
            # normalize image intensities
            tqdm_dask(mean := dask.persist(amps.mean(dim=['y','x'])), desc=f'Amplitude Normalization sw{subswath}')
            # dask.persist returns tuple
            norm = (mean[0]/mean[0].mean(dim='date'))
            del mean
            # compute Amplitude Dispersion Index (ADI)
            stats = (norm*amps).pipe(lambda x: (x.mean(dim='date'), x.std(dim='date')))
            del amps, norm
            ds = xr.merge([(scale*stats[0]).rename('amplitude'), (scale*stats[1]).rename('dispersion')])
            del stats
            self.save_model(ds.rename({'y': 'a', 'x': 'r'}), name=f'F{subswath}_adi', caption=f'Amplitude Dispersion Index (ADI) sw{subswath}')
            del ds

    def get_adi_threshold(self, subswath, threshold, chunksize=None):
        """
        Vectorize Amplitude Dispersion Index (ADI) raster values selected using the specified threshold.
        """
        import numpy as np
        import dask
        import pandas as pd
        import geopandas as gpd

        if chunksize is None:
            chunksize = self.chunksize
    
        def adi_block(ys, xs):
            from scipy.interpolate import griddata
            # we can calculate more accurate later
            dy = dx = 10
            trans_inv_block = trans_inv.sel(y=slice(min(ys)-dy,max(ys)+dy), x=slice(min(xs)-dx,max(xs)+dx))
            lt_block = trans_inv_block.lt.compute(n_workers=1).data.ravel()
            ll_block = trans_inv_block.ll.compute(n_workers=1).data.ravel()
            block_y, block_x = np.meshgrid(trans_inv_block.y.data, trans_inv_block.x.data)
            points = np.column_stack([block_y.ravel(), block_x.ravel()])
            # following NetCDF indices 0.5,1.5,...
            adi_block = adi.sel(y=slice(min(ys),max(ys)+1), x=slice(min(xs),max(xs)+1))
            adi_block_value = adi_block.compute(n_workers=1).data.ravel()
            adi_block_mask = adi_block_value<=threshold
            adi_block_value = adi_block_value[adi_block_mask]
            adi_block_y, adi_block_x = np.meshgrid(adi_block.y, adi_block.x)
            adi_block_y = adi_block_y.ravel()[adi_block_mask]
            adi_block_x = adi_block_x.ravel()[adi_block_mask]
            # interpolate geographic coordinates, coarsen=2 grid is required for the best accuracy
            grid_lt = griddata(points, lt_block, (adi_block_y, adi_block_x), method='linear').astype(np.float32)
            grid_ll = griddata(points, ll_block, (adi_block_y, adi_block_x), method='linear').astype(np.float32)
            # return geographic coordinates and values
            return np.column_stack([grid_lt, grid_ll, adi_block_value])
    
        # data grid and transform table
        adi = self.get_adi(subswath, chunksize=chunksize)
        trans_inv = self.get_trans_dat_inv(subswath, chunksize=chunksize)
    
        # split to equal chunks and rest
        ys_blocks = np.array_split(np.arange(adi.y.size), np.arange(0, adi.y.size, chunksize)[1:])
        xs_blocks = np.array_split(np.arange(adi.x.size), np.arange(0, adi.x.size, chunksize)[1:])
        # arrays size is unknown so we cannot construct dask array
        blocks = []
        for ys_block in ys_blocks:
            for xs_block in xs_blocks:
                block = dask.delayed(adi_block)(ys_block, xs_block)
                blocks.append(block)
                del block
    
        # materialize the result as a set of numpy arrays
        tqdm_dask(model := dask.persist(blocks), desc=f'Amplitude Dispersion Index (ADI) Threshold sw{subswath}')
        del blocks
        # the result is already calculated and compute() returns the result immediately
        model = np.concatenate(dask.compute(model)[0][0])
        # convert to geopandas object
        columns = {'adi': model[:,2], 'geometry': gpd.points_from_xy(model[:,1], model[:,0])}
        df = gpd.GeoDataFrame(columns, crs="EPSG:4326")
        del columns
        return df

    #df.to_file('adi.sqlite', driver='GPKG')
    def get_adi_threshold_parallel(self, threshold, **kwargs):
        subswaths = self.get_subswaths()
    
        adis = [self.get_adi_threshold(subswath, threshold, **kwargs) for subswath in subswaths]
        return adis[0] if len(adis)==1 else adis

    def geotif_stack(self, dates=None, intensity=False, chunksize=None):
        """
        tiffs = sbas.geotif_stack(['2022-06-16', '2022-06-28'], intensity=True)
        """
        import xarray as xr
        import rioxarray as rio
        import pandas as pd
        import numpy as np
        # from GMTSAR code
        DFACT = 2.5e-07
    
        if chunksize is None:
            chunksize = self.chunksize

        stack = []
        for subswath in self.get_subswaths():
            if dates is None:
                dates = self.df[self.df['subswath']==subswath].index.values
            #print ('dates', dates)
            tiffs = self.df[(self.df['subswath']==subswath)&(self.df.index.isin(dates))].datapath.values
            tiffs = [rio.open_rasterio(tif, chunks=chunksize)[0] for tif in tiffs]
            # build stack
            tiffs = xr.concat(tiffs, dim='date')
            tiffs['date'] = pd.to_datetime(dates)
            # 2 and 4 multipliers to have the same values as in SLC
            if intensity:
                stack.append((2*DFACT*np.abs(tiffs))**2)
            else:
                stack.append(2*DFACT*tiffs)

        return stack[0] if len(stack)==1 else stack

    def slc_stack(self, dates=None, intensity=False, chunksize=None):
        import xarray as xr
        import pandas as pd
        import numpy as np
        import dask

        if chunksize is None:
            chunksize = self.chunksize

        stack = []
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
            amps = [self.PRM(subswath, date).read_SLC_int(intensity=intensity, chunksize=minichunksize) for date in dates]
            # build stack
            amps = xr.concat(amps, dim='date')
            amps['date'] = pd.to_datetime(dates)
            stack.append(amps)

        return stack[0] if len(stack)==1 else stack

    def solid_tide(self, dates, coords, debug=False):
        """
        Compute the tidal correction using GMTSAR binary tool solid_tide.

        solid_tide yyyyddd.fffffff < lon_lat > lon_lat_dx_dy_dz

        coords = sbas.solid_tide(sbas.df.index, coords=[13.40076, 47.40143])
        coords = sbas.solid_tide(sbas.df.index, coords=[[13.40076, 47.40143]])
        coords = sbas.solid_tide(sbas.df.index, coords=[[13.40076, 47.40143], [13.40076, 47.40143]])

        """
        import numpy as np
        import pandas as pd
        from io import StringIO, BytesIO
        import os
        import subprocess
    
        coords = np.asarray(coords)
        if len(coords.shape) == 1:
            coords = [coords]
        buffer = BytesIO()
        np.savetxt(buffer, coords, delimiter=' ', fmt='%.6f')
        stdin_data = buffer.getvalue()
        #print ('stdin_data', stdin_data)

        outs = []
        for date in dates:
            SC_clock_start, SC_clock_stop = sbas.PRM(None, date).get('SC_clock_start', 'SC_clock_stop')
            dt = (SC_clock_start + SC_clock_stop)/2
            argv = ['solid_tide', str(dt)]
            #cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
            cwd = self.basedir
            p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, cwd=cwd, bufsize=10*1000*1000)
            stdout_data, stderr_data = p.communicate(input=stdin_data)

            stderr_data = stderr_data.decode('ascii')
            if stderr_data is not None and len(stderr_data) and debug:
                print ('DEBUG: solid_tide', stderr_data)
                return None

            out = np.fromstring(stdout_data, dtype=float, sep=' ')
            outs.append(out)
        return np.asarray(outs) if len(outs)>1 else outs[0]
