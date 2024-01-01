# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2022, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_landmask import Stack_landmask

class Stack_unwrap_snaphu(Stack_landmask):

    # -s for SMOOTH mode and -d for DEFO mode when DEFOMAX_CYCLE should be defined in the configuration
    # DEFO mode (-d) and DEFOMAX_CYCLE=0 is equal to SMOOTH mode (-s)
    # https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu_man1.html
    def snaphu(self, phase, corr=None, conf=None, conncomp=False, debug=False):
        """
        Unwraps phase using SNAPHU with the given phase and correlation data.

        This function unwraps the phase of an interferogram using the Statistical-cost, Network-flow Algorithm
        for Phase Unwrapping (SNAPHU) with user-defined parameters. The unwrapped phase is saved as a grid file
        in the working directory.

        Parameters
        ----------
        phase : xarray.DataArray
            The phase data as a string or xarray.DataArray, default is 'phasefilt'.

        corr : xarray.DataArray, optional
            The correlation data as a string or xarray.DataArray, default is 'corr'.

        conf : str, optional
            The SNAPHU configuration string, default is None (use the PRM's snaphu_config method).

        conncomp : bool, optional
            If True, return connection components map, default is False.

        debug : bool, optional
            If True, print debugging information during the unwrapping process, default is False.

        Returns
        -------
        xarray.Dataset
            Return the unwrapped phase and optional connection components as an xarray.Dataset.

        """
        import xarray as xr
        import numpy as np
        import pandas as pd
        import os
        import subprocess
        from datetime import datetime
        import uuid
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        from dask.distributed import utils_perf
        utils_perf.disable_gc_diagnosis()

        if conf is None:
            conf = self.snaphu_config()
        # set unique processing subdirectory
        conf += f'    TILEDIR snaphu_tiledir_{str(uuid.uuid4())}'

        # define basename for SNAPHU temp files
        # crop .grd from filename
        basename = self.get_filename(f'snaphu_{str(uuid.uuid4())}', '')[:-4]
        #print ('basename', basename)

        # SNAPHU input files
        phase_in = basename + '.phase'
        corr_in = basename + '.corr'
        mask_in = basename + '.bytemask'
        # SNAPHU output files
        unwrap_out = basename + 'unwrap.out'
        conncomp_out = basename + 'conncomp.out'

        # prepare SNAPHU input files
        # NaN values are not allowed for SNAPHU phase input file
        # interpolate when exist valid values around and fill zero pixels far away from valid ones
        # convert to lazy array setting the chunk size
        phase.fillna(0).compute(n_workers=1).values.astype(np.float32).tofile(phase_in)
        # SNAPHU masks out 0 and uses only valid pixels with mask 1
        xr.where(np.isnan(phase), 0, 1).values.astype(np.ubyte).tofile(mask_in)
    
        if corr is not None:
            # NaN values are not allowed for SNAPHU correlation input file
            # just fill NaNs by zeroes because the main trick is phase filling
            corr.fillna(0).compute(n_workers=1).values.astype(np.float32).tofile(corr_in)

        # launch SNAPHU binary (NaNs are not allowed for input but returned in output)
        argv = ['snaphu', phase_in, str(phase.shape[1]), '-M', mask_in, '-f', '/dev/stdin', '-o', unwrap_out, '-d']
        # output connection componetets map
        if conncomp:
            argv.append('-g')
            argv.append(conncomp_out)
        # add optional correlation grid
        if corr is not None:
            argv.append('-c')
            argv.append(corr_in)
        if debug:
            argv.append('-v')
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='utf8', bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=conf)

        outs = []
        # check for expected SNAPHU output files
        if os.path.exists(unwrap_out) and (not conncomp or os.path.exists(conncomp_out)):
            # convert to grid unwrapped phase from SNAPHU output applying postprocessing
            values = np.fromfile(unwrap_out, dtype=np.float32).reshape(phase.shape)
            #values = np.frombuffer(stdout_data, dtype=np.float32).reshape(phase.shape)
            # revert NaNs in output because SNAPNU does not support them
            unwrap = xr.DataArray(values, phase.coords, name='phase').chunk(self.chunksize).where(~np.isnan(phase))
            outs.append(unwrap)
            del values, unwrap

            if conncomp:
                # convert to grid the connected components from SNAPHU output as is (UCHAR)
                values = np.fromfile(conncomp_out, dtype=np.ubyte).reshape(phase.shape)
                conn = xr.DataArray(values, phase.coords, name='conncomp').chunk(self.chunksize)
                outs.append(conn)
                del values, conn
        else:
            # return the same data structure as expected but NaN-filled
            outs.append(xr.full_like(phase, np.nan).rename('phase'))
            if conncomp:
                outs.append(xr.full_like(phase, np.nan).rename('conncomp'))
        out = xr.merge(outs)
        del outs

        # add processing log
        out.attrs['snaphu'] = stdout_data + '\n' + stderr_data

        # the output files deleted immediately
        # but these are accessible while open descriptors persist
        for tmp_file in [phase_in, corr_in, mask_in, unwrap_out, conncomp_out]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

        return out

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
        #return self.PRM().snaphu_config(defomax, **kwargs)
        import os
        import joblib

        tiledir = os.path.splitext(self.PRM().filename)[0]
        n_jobs = joblib.cpu_count()

        conf_basic = f"""
        # basic config
        INFILEFORMAT   FLOAT_DATA
        OUTFILEFORMAT  FLOAT_DATA
        AMPFILEFORMAT  FLOAT_DATA
        CORRFILEFORMAT FLOAT_DATA
        ALTITUDE       693000.0
        EARTHRADIUS    6378000.0
        NEARRANGE      831000
        DR             18.4
        DA             28.2
        RANGERES       28
        AZRES          44
        LAMBDA         0.0554658
        NLOOKSRANGE    1
        NLOOKSAZ       1
        TILEDIR        {tiledir}_snaphu_tiledir
        NPROC          {n_jobs}
        """
        conf_custom = '# custom config\n'
        # defomax can be None
        keyvalues = ([('DEFOMAX_CYCLE', defomax)] if defomax is not None else []) + list(kwargs.items())
        for key, value in keyvalues:
            if isinstance(value, bool):
                value = 'TRUE' if value else 'FALSE'
            conf_custom += f'        {key} {value}\n'
        return conf_basic + conf_custom
