# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2022, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_landmask import SBAS_landmask

class SBAS_unwrap_snaphu(SBAS_landmask):

    # -s for SMOOTH mode and -d for DEFO mode when DEFOMAX_CYCLE should be defined in the configuration
    # DEFO mode (-d) and DEFOMAX_CYCLE=0 is equal to SMOOTH mode (-s)
    # https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu_man1.html
    def snaphu(self, phase, corr, conf=None, tiledir=None, chunksize=None, debug=False):
        """
        Unwraps phase using SNAPHU with the given phase and correlation data.

        This function unwraps the phase of an interferogram using the Statistical-cost, Network-flow Algorithm
        for Phase Unwrapping (SNAPHU) with user-defined parameters. The unwrapped phase is saved as a grid file
        in the working directory.

        Parameters
        ----------
        phase : str or xarray.DataArray, optional
            The phase data as a string or xarray.DataArray, default is 'phasefilt'.

        corr : str or xarray.DataArray, optional
            The correlation data as a string or xarray.DataArray, default is 'corr'.

        conf : str, optional
            The SNAPHU configuration string, default is None (use the PRM's snaphu_config method).

        interactive : bool, optional
            If True, return the unwrapped phase as an xarray.DataArray instead of saving it to disk, default is True.

        chunksize : tuple, optional
            The chunk size for dask arrays, default is None (use the instance's chunksize).

        debug : bool, optional
            If True, print debugging information during the unwrapping process, default is False.

        Returns
        -------
        xarray.DataArray
            If interactive is True, return the unwrapped phase as an xarray.DataArray; otherwise, return None and
            save the unwrapped phase to disk.

        """
        import xarray as xr
        import numpy as np
        import pandas as pd
        import os
        import subprocess

        # define lost class variables due to joblib
        if chunksize is None:
            chunksize = self.chunksize

        if conf is None:
            conf = self.PRM().snaphu_config()
        if tiledir is not None:
            conf += f'    TILEDIR {tiledir}'
        
        # define basename for SNAPHU temp files
        pair = phase.pair.replace(' ', '_')
        # crop .grd from filename
        basename = self.get_filenames(pair, '')[0][:-4]
        #print ('basename', basename)

        # SNAPHU input files
        phase_in = basename + 'unwrap.phase'
        corr_in = basename + 'unwrap.corr'
        # SNAPHU output files
        unwrap_out = basename + 'unwrap.out'
        conncomp_out = basename + 'conncomp.out'

        # cleanup from previous runs
        # conncomp_filename, unwrap_filename are required output files
        for tmp_file in [phase_in, corr_in, unwrap_out, conncomp_out]:
            #print ('tmp_file', tmp_file)
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

        # prepare SNAPHU input files
        # NaN values are not allowed for SNAPHU phase input file
        # interpolate when exist valid values around and fill zero pixels far away from valid ones
        self.nearest_grid(phase).fillna(0).astype(np.float32).values.tofile(phase_in)
        # NaN values are not allowed for SNAPHU correlation input file
        # just fill NaNs by zeroes because the main trick is phase filling
        corr.fillna(0).astype(np.float32).values.tofile(corr_in)

        # launch SNAPHU binary (NaNs are not allowed for input but returned in output)
        argv = ['snaphu', phase_in, str(phase.shape[1]), '-c', corr_in,
                '-f', '/dev/stdin', '-o', unwrap_out, '-g', conncomp_out, '-d']
        if debug:
            argv.append('-v')
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='ascii', bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=conf)

        # convert to grid the connected components from SNAPHU output as is (UCHAR)
        values = np.fromfile(conncomp_out, dtype=np.ubyte).reshape(phase.shape)
        conn = xr.DataArray(values, phase.coords, name='conncomp').chunk(chunksize)
        del values

        # convert to grid unwrapped phase from SNAPHU output applying postprocessing
        values = np.fromfile(unwrap_out, dtype=np.float32).reshape(phase.shape)
        #values = np.frombuffer(stdout_data, dtype=np.float32).reshape(phase.shape)
        unwrap = xr.DataArray(values, phase.coords, name='phase').chunk(chunksize)
        del values

        # the output files required for interactive output
        for tmp_file in [phase_in, corr_in] + [unwrap_out, conncomp_out] if not interactive else []:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

        if debug:
            return (unwrap, conn, stdout_data, stderr_data)
        else:
            # revert NaNs in output because SNAPNU does not support them
            return (unwrap.where(corr), conn)
            
    def snaphu_config(self, defomax=0, **kwargs):
        """
        Generate a Snaphu configuration file.

        Parameters
        ----------
        defomax : int, optional
            Maximum deformation value. Default is 0.
        **kwargs : dict, optional
            Additional parameters to include in the configuration file.

        Returns
        -------
        str
            The Snaphu configuration file content.

        Notes
        -----
        This method uses the `snaphu_config` method of the PRM object.

        Examples
        --------
        Generate a Snaphu configuration file with defomax=10:
        snaphu_config(defomax=10)

        Generate a Snaphu configuration file with defomax=5 and additional parameters:
        snaphu_config(defomax=5, param1=10, param2=20)
        """
        return self.PRM().snaphu_config(defomax, **kwargs)
 