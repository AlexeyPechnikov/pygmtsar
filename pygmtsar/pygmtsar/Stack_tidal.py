# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_incidence import Stack_incidence

class Stack_tidal(Stack_incidence):

#     # naive code is not suitable for large grids and stacks
#     def cube_tidal_los(self, data):
#         """
#         Calculate tidal LOS displacement [m] for data dates and spatial extent
#         """
#         import xarray as xr
# 
#         # suppose only one subswath
#         solid_tide = self.get_tidal()[0]
#         # interpolate on the data_pairs 2D grid
#         tidal_dates = solid_tide.interp_like(data, method='linear', assume_sorted=True)
#         # compute LOS projection, [m]
#         tidal_dates_los = self.los_projection(tidal_dates)
#         return tidal_dates_los

    def cube_tidal_los(self, data_pairs, chunksize=None):
        """
        Calculate tidal LOS displacement [m] for data pairs dates and spatial extent
        """
        import pandas as pd
        import dask
        import xarray as xr
        import numpy as np

        if chunksize is None:
            chunksize = self.chunksize

        # select only required dates
        dates = self.get_pairs(data_pairs, dates=True)[1]
        data_grid = data_pairs[0] if len(data_pairs.dims) == 3 else data_pairs
    
        # suppose only one subswath
        solid_tide = self.get_tidal()[0].sel(date=dates)

        # satellite look vector
        sat_look = self.get_sat_look()

        @dask.delayed
        def interp_block(date, ys_block, xs_block):
            # use outer variables
            # interpolate on the data_pairs 2D grid
            coords = {'y': ys_block, 'x': xs_block}
            block_tidal = solid_tide.sel(date=date).interp(coords, method='linear', assume_sorted=True)
            block_look = sat_look.interp(coords, method='linear', assume_sorted=True)
            los = xr.dot(xr.concat([block_look.look_E, block_look.look_N, block_look.look_U], dim='dim'),
                      xr.concat([block_tidal.dx, block_tidal.dy, block_tidal.dz], dim='dim'),
                      dims=['dim'])
            del block_tidal, block_look, coords
            return los.data[None,]

        # define output radar coordinates grid and split to equal chunks and rest
        ys_blocks = np.array_split(data_pairs.y, np.arange(0, data_pairs.y.size, chunksize)[1:])
        xs_blocks = np.array_split(data_pairs.x, np.arange(0, data_pairs.x.size, chunksize)[1:])

        # per-block processing
        blocks3d  = []
        for date in dates:
            blocks2d  = []
            for ys_block in ys_blocks:
                blocks = []
                for xs_block in xs_blocks:
                    block = dask.array.from_delayed(interp_block(date, ys_block, xs_block),
                                                    shape=(1, ys_block.size, xs_block.size), dtype=np.float32)
                    blocks.append(block)
                    del block
                blocks2d.append(blocks)
                del blocks
            blocks3d.append(blocks2d)
            del blocks2d
        dask_block = dask.array.block(blocks3d)
        del blocks3d

        # convert coordinate to valid dates
        coords = {'date': pd.to_datetime(dates), 'y': data_grid.coords['y'], 'x': data_grid.coords['x']}
        out = xr.DataArray(dask_block, coords=coords)
        del dask_block
        return out.rename('tidal')

    def cube_tidal_los_rad(self, data, chunksize=None):
        """
        Calculate tidal LOS displacement [rad] for data dates and spatial extent
        """
        return 1000*self.cube_tidal_los(data, chunksize=chunksize)/self.los_displacement_mm(1)

    def stack_tidal_los(self, data_pairs, chunksize=None):
        """
        Calculate tidal LOS displacement [m] for data_pairs pairs and spatial extent
        """
        import xarray as xr

        # extract pairs
        pairs = self.get_pairs(data_pairs)
        # interpolate on the data_pairs 2D grid
        # compute LOS projection, [m]
        tidal_dates_los = self.cube_tidal_los(data_pairs, chunksize=chunksize)
        # calculate differences between end and start dates for all the pairs
        tidal_pairs = []
        for rec in pairs.itertuples():
            tidal_pair = tidal_dates_los.sel(date=rec.rep) - tidal_dates_los.sel(date=rec.ref)
            tidal_pairs.append(tidal_pair)
        # form 3D stack
        return xr.concat(tidal_pairs, dim='pair').assign_coords({'pair': data_pairs.pair})

#     # faster, but there is no significant advantage vs naive version
#     def stack_tidal_los(self, data_pairs, chunksize=None):
#         """
#         Calculate tidal LOS displacement [m] for data pairs dates and spatial extent
#         """
#         import pandas as pd
#         import dask
#         import xarray as xr
#         import numpy as np
# 
#         if chunksize is None:
#             chunksize = self.chunksize
# 
#         # select only required dates
#         pairs, dates = self.get_pairs(data_pairs, dates=True)
#         data_grid = data_pairs[0] if len(data_pairs.dims) == 3 else data_pairs
# 
#         # suppose only one subswath
#         solid_tide = self.get_tidal()[0].sel(date=dates)
# 
#         # satellite look vector
#         sat_look = self.get_sat_look()
# 
#         @dask.delayed
#         def interp_block(date_ref, date_rep, ys_block, xs_block):
#             # use outer variables
#             # interpolate on the data_pairs 2D grid
#             coords = {'y': ys_block, 'x': xs_block}
#             block_tidal = solid_tide.sel(date=date_rep).interp(coords, method='linear', assume_sorted=True) - \
#                           solid_tide.sel(date=date_ref).interp(coords, method='linear', assume_sorted=True)
#             block_look = sat_look.interp(coords, method='linear', assume_sorted=True)
#             los = xr.dot(xr.concat([block_look.look_E, block_look.look_N, block_look.look_U], dim='dim'),
#                       xr.concat([block_tidal.dx, block_tidal.dy, block_tidal.dz], dim='dim'),
#                       dims=['dim'])
#             del block_tidal, block_look, coords
#             return los.data[None,]
# 
#         # define output radar coordinates grid and split to equal chunks and rest
#         ys_blocks = np.array_split(data_pairs.y, np.arange(0, data_pairs.y.size, chunksize)[1:])
#         xs_blocks = np.array_split(data_pairs.x, np.arange(0, data_pairs.x.size, chunksize)[1:])
# 
#         # per-block processing
#         blocks3d  = []
#         for rec in pairs.itertuples():
#             blocks2d  = []
#             for ys_block in ys_blocks:
#                 blocks = []
#                 for xs_block in xs_blocks:
#                     block = dask.array.from_delayed(interp_block(rec.ref, rec.rep, ys_block, xs_block),
#                                                     shape=(1, ys_block.size, xs_block.size), dtype=np.float32)
#                     blocks.append(block)
#                     del block
#                 blocks2d.append(blocks)
#                 del blocks
#             blocks3d.append(blocks2d)
#             del blocks2d
#         dask_block = dask.array.block(blocks3d)
#         del blocks3d
# 
#         out = xr.DataArray(dask_block, coords=data_pairs.coords)
#         del dask_block
#         return out.rename('tidal')

    def stack_tidal_los_rad(self, data_pairs, chunksize=None):
        """
        Calculate tidal LOS displacement [rad] for data_pairs pairs and spatial extent
        """
        return 1000*self.stack_tidal_los(data_pairs, chunksize=chunksize)/self.los_displacement_mm(1)

    def stack_tidal_correction_wrap(self, data_pairs, chunksize=None):
        """
        Apply tidal correction to wrapped phase pairs [rad] and wrap the result.
        """
        return self.wrap(data_pairs - self.stack_tidal_los_rad(data_pairs, chunksize=chunksize))
    
    def get_tidal(self, subswath=None, chunksize=None):
        return self.open_grid('tidal', subswath=subswath, chunksize=chunksize)

    def tidal_parallel(self, pairs, coarsen=32, chunksize=None, interactive=False):
        import xarray as xr
        import pandas as pd
        import numpy as np
        from io import StringIO, BytesIO
        import subprocess
        #import os

        # some scenes can be missed in pairs 
        dates = self.get_pairs(pairs, dates=True)[1]

        # expand simplified definition
        if not isinstance(coarsen, (list,tuple, np.ndarray)):
            coarsen = (coarsen, coarsen)

        subswath = self.get_subswath()
        trans_inv = self.get_trans_inv(subswath)
        dy, dx = np.diff(trans_inv.y)[0], np.diff(trans_inv.x)[0]
        step_y, step_x = int(np.round(coarsen[0]*dy)), int(np.round(coarsen[1]*dx))
        #print ('step_y, step_x', step_y, step_x)
        grid = trans_inv.sel(y=trans_inv.y[step_y//2::step_y], x=trans_inv.x[step_x//2::step_x])

        coords = np.column_stack([grid.ll.values.ravel(), grid.lt.values.ravel()])    
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
            if stderr_data is not None and len(stderr_data):
                print ('DEBUG: solid_tide', stderr_data)
                return None

            out = np.fromstring(stdout_data, dtype=float, sep=' ')
            outs.append(out.reshape(grid.y.size, grid.x.size, 5))

        outs = np.asarray(outs)
        coords = {'date': pd.to_datetime(dates), 'y': grid.y, 'x': grid.x}
        das = {v: xr.DataArray(outs[...,idx], coords=coords) for (idx, v) in enumerate(['lon', 'lat', 'dx', 'dy', 'dz'])}
        ds = xr.Dataset(das)

        if interactive:
                return ds
        #.rename({'y': 'a', 'x': 'r'})
        return self.save_grid(ds, 'tidal', subswath, f'Solid Earth Tides Computing sw{subswath}', chunksize)

#     def solid_tide(self, dates, data, debug=False):
#         """
#         Compute the tidal correction using the GMTSAR binary tool solid_tide for geographic coordinate points.
# 
#         Parameters
#         ----------
#         dates : array_like
#             Dates for which to compute the tidal correction. Should be in the format yyyyddd.fffffff.
#         coords : array_like
#             Coordinates for which to compute the tidal correction. Can be a single pair of coordinates [lon, lat],
#             or a list of pairs [[lon1, lat1], [lon2, lat2], ...].
#         debug : bool, optional
#             If True, print debug information. Default is False.
# 
#         Returns
#         -------
#         pandas.DataFrame
#             DataFrame with columns 'lon', 'lat', 'dx', 'dy', 'dz', indexed by 'date'.
# 
#         Examples
#         --------
#         Compute the tidal correction for a single pair of coordinates:
# 
#             coords = stack.solid_tide(stack.df.index, coords=[13.40076, 47.40143])
# 
#         Compute the tidal correction for multiple pairs of coordinates:
# 
#             coords = stack.solid_tide(stack.df.index, coords=[[13.40076, 47.40143], [13.40076, 47.40143]])
#     
#         Compute the tidal correction for point geodataframe:
#             coords = stack.solid_tide(stack.df.index, AOI)
# 
#         Compute the tidal correction for a single record point geodataframe:
#             coords = stack.solid_tide(stack.df.index, AOI.head(1))
# 
#         Output:
# 
#             >>> stack.solid_tide(stack.df.index[:3], coords=[lon, lat])        
#             lon	lat	dx	dy	dz
#             date					
#             2022-06-16	13.400758	47.401431	-0.066918	-0.004765	0.016200
#             2022-06-28	13.400758	47.401431	-0.033571	0.012279	-0.099899
#             2022-07-10	13.400758	47.401431	-0.000806	-0.007983	-0.150675
# 
#             >>> stack.solid_tide(stack.df.index[:3], coords=[[lon, lat], [lon+1, lat+1]])
#             lon	lat	dx	dy	dz
#             date					
#             2022-06-16	13.400758	47.401431	-0.066918	-0.004765	0.016200
#             2022-06-16	14.400758	48.401431	-0.066594	-0.004782	0.010346
#             2022-06-28	13.400758	47.401431	-0.033571	0.012279	-0.099899
#             2022-06-28	14.400758	48.401431	-0.033011	0.012323	-0.100986
#             2022-07-10	13.400758	47.401431	-0.000806	-0.007983	-0.150675
#             2022-07-10	14.400758	48.401431	-0.001107	-0.007758	-0.151605
# 
#         Notes
#         -----
#         This function computes the tidal correction in geographic coordinates based on the dates and coordinates provided.
#         The correction is computed by calling the GMTSAR binary tool 'solid_tide' with the date and coordinates as input.
#         """
#         import numpy as np
#         import geopandas as gpd
#         import pandas as pd
#         from io import StringIO, BytesIO
#         #import os
#         import subprocess
# 
#         if isinstance(data, gpd.GeoDataFrame):
#             coords = np.array([(geom.x, geom.y) for geom in data.geometry])
#         else:
#             coords = np.asarray(data)
#             if len(coords.shape) == 1:
#                 coords = [coords]
#         #print ('coords', coords)
#         buffer = BytesIO()
#         np.savetxt(buffer, coords, delimiter=' ', fmt='%.6f')
#         stdin_data = buffer.getvalue()
#         #print ('stdin_data', stdin_data)
# 
#         outs = []
#         for date in dates:
#             SC_clock_start, SC_clock_stop = self.PRM(None, date).get('SC_clock_start', 'SC_clock_stop')
#             dt = (SC_clock_start + SC_clock_stop)/2
#             argv = ['solid_tide', str(dt)]
#             #cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
#             cwd = self.basedir
#             p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
#                                  stderr=subprocess.PIPE, cwd=cwd, bufsize=10*1000*1000)
#             stdout_data, stderr_data = p.communicate(input=stdin_data)
# 
#             stderr_data = stderr_data.decode('ascii')
#             if stderr_data is not None and len(stderr_data) and debug:
#                 print ('DEBUG: solid_tide', stderr_data)
#                 return None
# 
#             out = np.fromstring(stdout_data, dtype=float, sep=' ')
#             outs.append(out)
# 
#         return pd.DataFrame(np.asarray(outs).reshape(-1,5),
#                             columns=['lon', 'lat', 'dx', 'dy', 'dz'],
#                             index=np.repeat(dates, len(coords)))

    # use this direct calculation function to check the new version based on trans_inv (incorrect for now)
    def solid_tide2(self, dates, data, debug=False):
        """
        Compute the tidal correction using the GMTSAR binary tool solid_tide for geographic coordinate points.

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

            coords = stack.solid_tide(stack.df.index, coords=[13.40076, 47.40143])

        Compute the tidal correction for multiple pairs of coordinates:

            coords = stack.solid_tide(stack.df.index, coords=[[13.40076, 47.40143], [13.40076, 47.40143]])

        Compute the tidal correction for point geodataframe:
            coords = stack.solid_tide(stack.df.index, AOI)

        Compute the tidal correction for a single record point geodataframe:
            coords = stack.solid_tide(stack.df.index, AOI.head(1))

        Output:

            >>> stack.solid_tide(stack.df.index[:3], coords=[lon, lat])        
            lon	lat	dx	dy	dz
            date					
            2022-06-16	13.400758	47.401431	-0.066918	-0.004765	0.016200
            2022-06-28	13.400758	47.401431	-0.033571	0.012279	-0.099899
            2022-07-10	13.400758	47.401431	-0.000806	-0.007983	-0.150675

            >>> stack.solid_tide(stack.df.index[:3], coords=[[lon, lat], [lon+1, lat+1]])
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
        import geopandas as gpd
        import pandas as pd
        from io import StringIO, BytesIO
        #import os
        import subprocess

        if isinstance(data, gpd.GeoDataFrame):
            coords = np.array([(geom.x, geom.y) for geom in data.geometry])
        else:
            coords = np.asarray(data)
            if len(coords.shape) == 1:
                coords = [coords]
        #print ('coords', coords)
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
