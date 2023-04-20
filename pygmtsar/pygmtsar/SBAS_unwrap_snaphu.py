#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar
from .SBAS_landmask import SBAS_landmask

class SBAS_unwrap_snaphu(SBAS_landmask):

    # -s for SMOOTH mode and -d for DEFO mode when DEFOMAX_CYCLE should be defined in the configuration
    # DEFO mode (-d) and DEFOMAX_CYCLE=0 is equal to SMOOTH mode (-s)
    # https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu_man1.html
    def unwrap(self, pair, threshold=1e-3, conf=None, func=None, mask=None, conncomp=False,
               phase='phasefilt', corr='corr', interactive=True, chunksize=None, debug=False):
        import xarray as xr
        import numpy as np
        import os
        import subprocess

        # define lost class variables due to joblib
        if chunksize is None:
            chunksize = self.chunksize

        if conf is None:
            conf = self.PRM().snaphu_config()
    
        # open input data grids if needed
        if isinstance(phase, str):
            phase = self.open_grids([pair], phase, chunksize=chunksize, interactive=False)[0]
        if isinstance(corr, str):
            corr = self.open_grids([pair], corr, chunksize=chunksize, interactive=False)[0]
        if mask is not None and isinstance(mask, str):
            mask = self.open_grids(None, mask, chunksize=chunksize, interactive=False)
        
        # convert user-defined mask to boolean mask (NaN values converted to 0)
        if mask is not None:
            assert self.is_ra(mask), 'ERROR: mask should be defined in radar coordinates'
            # define mask on the same grid as phase and convert to boolean
            binmask = xr.where(mask>0, 1, 0).astype(bool)
            # unify grids to maybe smallest mask
            phase = phase.reindex_like(mask)
            corr = corr.reindex_like(mask)
        else:
            binmask = 1
        # add threshold mask
        binmask = binmask & (corr>=threshold)

        # extract dates from pair
        date1, date2 = pair

        basename = self.get_filenames(None, [pair], '')[0][:-4]
        #print ('basename', basename)

        # output data grids
        unwrap_filename = basename + 'unwrap.grd'
        conncomp_filename = basename + 'conncomp.grd'
        # SNAPHU input files
        phase_in = basename + 'unwrap.phase'
        corr_in = basename + 'unwrap.corr'
        # SNAPHU output files
        unwrap_out = basename + 'unwrap.out'
        conncomp_out = basename + 'conncomp.out'

        # cleanup from previous runs
        for tmp_file in [phase_in, corr_in, unwrap_out, conncomp_out] \
                      + [conncomp_filename, unwrap_filename] if not interactive else []:
            #print ('tmp_file', tmp_file)
            if os.path.exists(tmp_file):
                if debug:
                    print ('DEBUG: remove', tmp_file)
                os.remove(tmp_file)

        # prepare SNAPHU input files
        # NaN values are not allowed for SNAPHU phase input file
        # interpolate when exist valid values around and fill zero pixels far away from valid ones
        self.nearest_grid(phase.where(binmask)).fillna(0).astype(np.float32).values.tofile(phase_in)
        # NaN values are not allowed for SNAPHU correlation input file
        # just fill NaNs by zeroes because the main trick is phase filling
        corr.where(binmask).fillna(0).astype(np.float32).values.tofile(corr_in)

        # launch SNAPHU binary (NaNs are not allowed for input but returned in output)
        argv = ['snaphu', phase_in, str(phase.shape[1]), '-c', corr_in,
                '-f', '/dev/stdin', '-o', unwrap_out, '-d']
        if conncomp:
            argv.append('-g')
            argv.append(conncomp_out)
        if debug:
            argv.append('-v')
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='ascii', bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=conf)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: snaphu', stderr_data)
        if debug:
            print ('DEBUG: snaphu', stdout_data)

        if conncomp:
            # convert to grid the connected components from SNAPHU output as is (UCHAR)
            values = np.fromfile(conncomp_out, dtype=np.ubyte).reshape(phase.shape)
            conn = xr.DataArray(values, phase.coords, name='conncomp').chunk(chunksize)

        # convert to grid unwrapped phase from SNAPHU output applying postprocessing
        values = np.fromfile(unwrap_out, dtype=np.float32).reshape(phase.shape)
        #values = np.frombuffer(stdout_data, dtype=np.float32).reshape(phase.shape)
        unwrap = xr.DataArray(values, phase.coords).chunk(chunksize)
        # apply user-defined function for post-processing
        if func is not None:
            # the both grids should have the right chunksize
            unwrap = func(corr.where(binmask), unwrap)
        # apply binary mask after the post-processing to completely exclude masked regions
        # NaN values allowed in the output grid, assign userfriendly name for the output grid
        unwrap = unwrap.where(binmask).rename('phase')

        # the output files required for interactive output
        for tmp_file in [phase_in, corr_in] + [unwrap_out, conncomp_out] if not interactive else []:
            if os.path.exists(tmp_file):
                if debug:
                    print ('DEBUG: remove', tmp_file)
                os.remove(tmp_file)

        if interactive:
            return (unwrap, conn) if conncomp else unwrap

        # not interactive mode, save all the results to disk
        if conncomp:
            conn.to_netcdf(conncomp_filename, encoding={'conncomp': self.compression(chunksize=chunksize)}, engine=self.engine)
        # save to NetCDF file
        unwrap.to_netcdf(unwrap_filename, encoding={'phase': self.compression(chunksize=chunksize)}, engine=self.engine)
