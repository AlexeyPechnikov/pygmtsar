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

    def get_ps(self, subswath=None, chunksize=None):
        return self.open_grid('ps', subswath, chunksize=chunksize)

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
    def ps_parallel(self, dates=None, intensity=True, dfact=2.5e-07, chunksize=None):
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
            # intensity=False means complex data
            slcs = self.open_stack_slc(dates=dates, subswath=subswath, intensity=True, dfact=dfact)
            # convert to amplitude for GMTSAR compatible calculations
            if not intensity:
                slcs = np.sqrt(slcs)
            # normalize image intensities
            tqdm_dask(mean := dask.persist(slcs.mean(dim=['y','x'])), desc=f'Amplitude Normalization sw{subswath}')
            # dask.persist returns tuple
            norm = mean[0].mean(dim='date') / mean[0]
            # compute average and std.dev.
            stats = (norm * slcs).pipe(lambda x: (x.mean(dim='date'), x.std(dim='date')))
            del slcs
            ds = xr.merge([stats[0].rename('average'), stats[1].rename('deviation'), mean[0].rename('stack_average')])
            del stats, norm
            self.save_grid(ds.rename({'y': 'a', 'x': 'r'}), 'ps', subswath, f'Persistent Scatterers sw{subswath}')
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

    def solid_tide(self, dates, coords, debug=False):
        """
        Compute the tidal correction using the GMTSAR binary tool solid_tide.

        Parameters
        ----------
        dates : array_like
            Dates for which to compute the tidal correction. Should be in the format yyyyddd.fffffff.
        coords : array_like
            Coordinates for which to compute the tidal correction. Can be a single pair of coordinates [lon, lat],
            or a list of pairs [[lon1, lat1], [lon2, lat2], ...].
        debug : bool, optional
            If True, print debug information. Default is False.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'lon', 'lat', 'dx', 'dy', 'dz', indexed by 'date'.

        Examples
        --------
        Compute the tidal correction for a single pair of coordinates:

            coords = sbas.solid_tide(sbas.df.index, coords=[13.40076, 47.40143])

        Compute the tidal correction for multiple pairs of coordinates:

            coords = sbas.solid_tide(sbas.df.index, coords=[[13.40076, 47.40143], [13.40076, 47.40143]])

        Output:

            >>> sbas.solid_tide(sbas.df.index[:3], coords=[lon, lat])        
            lon	lat	dx	dy	dz
            date					
            2022-06-16	13.400758	47.401431	-0.066918	-0.004765	0.016200
            2022-06-28	13.400758	47.401431	-0.033571	0.012279	-0.099899
            2022-07-10	13.400758	47.401431	-0.000806	-0.007983	-0.150675

            >>> sbas.solid_tide(sbas.df.index[:3], coords=[[lon, lat], [lon+1, lat+1]])
            lon	lat	dx	dy	dz
            date					
            2022-06-16	13.400758	47.401431	-0.066918	-0.004765	0.016200
            2022-06-16	14.400758	48.401431	-0.066594	-0.004782	0.010346
            2022-06-28	13.400758	47.401431	-0.033571	0.012279	-0.099899
            2022-06-28	14.400758	48.401431	-0.033011	0.012323	-0.100986
            2022-07-10	13.400758	47.401431	-0.000806	-0.007983	-0.150675
            2022-07-10	14.400758	48.401431	-0.001107	-0.007758	-0.151605

        Notes
        -----
        This function computes the tidal correction in geographic coordinates based on the dates and coordinates provided.
        The correction is computed by calling the GMTSAR binary tool 'solid_tide' with the date and coordinates as input.
        """
        import numpy as np
        import pandas as pd
        from io import StringIO, BytesIO
        #import os
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
            SC_clock_start, SC_clock_stop = self.PRM(None, date).get('SC_clock_start', 'SC_clock_stop')
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

        return pd.DataFrame(np.asarray(outs).reshape(-1,5),
                            columns=['lon', 'lat', 'dx', 'dy', 'dz'],
                            index=np.repeat(dates, len(coords)))