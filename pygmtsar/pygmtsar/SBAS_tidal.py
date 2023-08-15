# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_incidence import SBAS_incidence

class SBAS_tidal(SBAS_incidence):

    def solid_tide(self, dates, data, debug=False):
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

            coords = sbas.solid_tide(sbas.df.index, coords=[13.40076, 47.40143])

        Compute the tidal correction for multiple pairs of coordinates:

            coords = sbas.solid_tide(sbas.df.index, coords=[[13.40076, 47.40143], [13.40076, 47.40143]])
    
        Compute the tidal correction for point geodataframe:
            coords = sbas.solid_tide(sbas.df.index, AOI)

        Compute the tidal correction for a single record point geodataframe:
            coords = sbas.solid_tide(sbas.df.index, AOI.head(1))

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
