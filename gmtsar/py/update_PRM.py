#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to read, write and update PRM files and calculate Doppler orbit by calc_dop_orb command line tool
#import pytest

class PRM:

    # replacement function for GMT based robust 2D trend coefficient calculations:
    # gmt trend2d r.xyz -Fxyzmw -N1r -V
    # gmt trend2d r.xyz -Fxyzmw -N2r -V
    # gmt trend2d r.xyz -Fxyzmw -N3r -V
    # https://github.com/GenericMappingTools/gmt/blob/master/src/trend2d.c#L719-L744
    # 3 model parameters
    # rank = 3 => nu = size-3
    @staticmethod
    def GMT_trend2d(data, rank):
        import numpy as np
        from sklearn.linear_model import LinearRegression
        # scale factor for normally distributed data is 1.4826
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.median_abs_deviation.html
        MAD_NORMALIZE = 1.4826
        # significance value
        sig_threshold = 0.51

        if rank not in [1,2,3]:
            raise Exception('Number of model parameters "rank" should be 1, 2, or 3')

        #see gmt_stat.c
        def gmtstat_f_q (chisq1, nu1, chisq2, nu2):
            import scipy.special as sc

            if chisq1 == 0.0:
                return 1
            if chisq2 == 0.0:
                return 0
            return sc.betainc(0.5*nu2, 0.5*nu1, chisq2/(chisq2+chisq1))

        if rank in [2,3]:
            x = data[:,0]
            x = np.interp(x, (x.min(), x.max()), (-1, +1))
        if rank == 3:
            y = data[:,1]
            y = np.interp(y, (y.min(), y.max()), (-1, +1))
        z = data[:,2]
        w = np.ones(z.shape)

        if rank == 1:
            xy = np.expand_dims(np.zeros(z.shape),1)
        elif rank == 2:
            xy = np.expand_dims(x,1)
        elif rank == 3:
            xy = np.stack([x,y]).transpose()

        # create linear regression object
        mlr = LinearRegression()

        chisqs = []
        coeffs = []
        while True:
            # fit linear regression
            mlr.fit(xy, z, sample_weight=w)

            r = np.abs(z - mlr.predict(xy))
            chisq = np.sum((r**2*w))/(z.size-3)    
            chisqs.append(chisq)
            k = 1.5 * MAD_NORMALIZE * np.median(r)
            w = np.where(r <= k, 1, (2*k/r) - (k * k/(r**2)))
            sig = 1 if len(chisqs)==1 else gmtstat_f_q(chisqs[-1], z.size-3, chisqs[-2], z.size-3)
            # Go back to previous model only if previous chisq < current chisq
            if len(chisqs)==1 or chisqs[-2] > chisqs[-1]:
                coeffs = [mlr.intercept_, *mlr.coef_]

            #print ('chisq', chisq, 'significant', sig)
            if sig < sig_threshold:
                break

        # get the slope and intercept of the line best fit
        return (coeffs[:rank])

    @staticmethod
    def from_list(prm_list):
        from io import StringIO
        prm = StringIO('\n'.join(prm_list))
        return PRM._from_io(prm)

    @staticmethod
    def from_str(prm_string):
        from io import StringIO
        prm = StringIO(prm_string)
        return PRM._from_io(prm)

    @staticmethod
    def from_file(prm_filename):
        #data = json.loads(document)
        prm = PRM._from_io(prm_filename)
        prm.filename = prm_filename
        return prm

    @staticmethod
    def _from_io(prm):
        import pandas as pd
        return PRM(pd.read_csv(prm, sep='\s+=\s+', header=None, names=['name', 'value'], engine='python').set_index('name')\
                    .applymap(lambda val : pd.to_numeric(val,errors='ignore')))

    def __init__(self, prm=None):
        import pandas as pd
        #print ('__init__')
        if prm is None:
            _prm = pd.DataFrame(None,columns=['name','value'])
        elif isinstance(prm, pd.DataFrame):
            _prm = prm.reset_index()
        else:
            _prm = prm.df.reset_index()
        self.df = _prm[['name', 'value']].drop_duplicates(keep='last', inplace=False).set_index('name')\
            .applymap(lambda val : pd.to_numeric(val,errors='ignore'))
        self.filename = None

    def __eq__(self, other):
        return isinstance(self, PRM) and self.df == other.df

    def __str__(self):
        return self.to_str()

    def __repr__(self):
        if self.filename:
            return 'Object %s (%s) %d items\n%r' % (self.__class__.__name__, self.filename, len(self.df), self.df)
        else:
            return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    # use 'g' format for Python and numpy float values
    def set(self, prm=None, gformat=False, **kwargs):
        import numpy as np

        if isinstance(prm, PRM):
            for (key, value) in prm.df.itertuples():
                self.df.loc[key] = value
        elif prm is not None:
            raise Exception('Arguments is not a PRM object')
        for key, value in kwargs.items():
            self.df.loc[key] = float(format(value, 'g')) if gformat and type(value) \
                in [float, np.float16, np.float32, np.float64, np.float128] else value
        return self

    def to_dataframe(self):
        return self.df

    def to_file(self, prm):
        self._to_io(prm)
        return self

    #def update(self):
    #    if self.filename is None:
    #        raise Exception('PRM is not created from file, use to_file() method instead')
    #    return self._to_io(self.filename)
    def update(self, name=None, safe=False, debug=False):
        """
        Save PRM file to disk and rename all 3 linked files together: input, LED, SLC if "name" defined
        If safe=True, save old (small) PRM and LED files and move only (large) SLC file.
        """
        import shutil
        import os

        if self.filename is None:
            raise Exception('PRM is not created from file, use to_file() method instead')

        if debug:
            print ('Debug mode: only print expected operations but do not perform the actual job')

        # rename linked files
        if name is not None:
            # define current files directory
            dirname0 = os.path.dirname(self.filename)
            # define new files basename
            basename = os.path.splitext(name)[0]
            shortname = os.path.split(basename)[-1]

            # rename PRM file
            if not safe:
                print (f'Remove old PRM file {self.filename} and save new one {name}')
                if not debug:
                    os.remove(self.filename)
            else:
                print (f'Safe mode: remain old PRM file {self.filename} and save new one {name}')
            # will be saved later
            self.filename = name

            input_file0 = os.path.join(dirname0, os.path.split(self.get('input_file'))[-1])
            input_ext = os.path.splitext(input_file0)[1]
            input_file = basename + input_ext
            print (f'PRM change input_file attribute {input_file0} -> {input_file}')
            self.set(input_file = f'{shortname}{input_ext}')
            if os.path.isfile(input_file0) and not input_file == input_file0:
                if not safe:
                    print (f'Rename input_file {input_file0} -> {input_file}')
                    if not debug:
                        os.rename(input_file0, input_file)
                else:
                    print (f'Safe mode: copy input_file {input_file0} -> {input_file}')
                    if not debug:
                        shutil.copy2(input_file0, input_file, follow_symlinks=True)

            SLC_file0 = os.path.join(dirname0, os.path.split(self.get('SLC_file'))[-1])
            SLC_file = basename + '.SLC'
            #print ('SLC_file', SLC_file)
            print (f'PRM change SLC_file attribute {SLC_file0} -> {SLC_file}')
            self.set(SLC_file = f'{shortname}.SLC')
            if os.path.isfile(SLC_file0) and not SLC_file == SLC_file0:
                print (f'Rename SLC_file {SLC_file0} -> {SLC_file}')
                if not debug:
                    os.rename(SLC_file0, SLC_file)

            led_file0 = os.path.join(dirname0, os.path.split(self.get('led_file'))[-1])
            led_file = basename + '.LED'
            #print ('led_file', led_file)
            print (f'PRM change LED_file attribute {led_file0} -> {led_file}')
            self.set(led_file = f'{shortname}.LED')
            if os.path.isfile(led_file0) and not led_file == led_file0:
                if not safe:
                    print (f'Rename LED_file {led_file0} -> {led_file}')
                    if not debug:
                        os.rename(led_file0, led_file)
                else:
                    print (f'Safe mode: copy LED_file {led_file0} -> {led_file}')
                    if not debug:
                        shutil.copy2(led_file0, led_file)

        if debug:
            return self
            #.sel('input_file','SLC_file','led_file')

        return self.to_file(self.filename)

    def to_str(self):
        return self._to_io()

    def _to_io(self, output=None):
        return self.df.reset_index().astype(str).apply(lambda row: (' = ').join(row), axis=1)\
            .to_csv(output, header=None, index=None)

    def sel(self, *args):
        return PRM(self.df.loc[[*args]])

    def __add__(self, other):
        import pandas as pd
        if isinstance(other, PRM):
            prm = pd.concat([self.df, other.df])
            # drop duplicates
            prm = prm.groupby(prm.index).last()
        else:
            prm = self.df + other
        return PRM(prm)

    def get(self, *args):
        out = [self.df.loc[[key]].iloc[0].values[0] for key in args]
        if len(out) == 1:
            return out[0]
        return out

    def calc_dop_orb(self, earth_radius=0, doppler_centroid=0, inplace=False):
        """
        Usage: calc_dop_orb  file.PRM  added.PRM  earth_radius  [doppler_centroid]
            file.PRM     - input name of PRM file 
            new.PRM      - output additional parameters to add to the PRM file 
            earth_radius - input set earth radius, 0 calculates radius 
            [doppler_centroid] - no parameter calculates doppler 
            [doppler_centroid] - use value (e.g. 0.0) to force doppler 
        """
        import subprocess
        import os
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(['calc_dop_orb', '/dev/stdin', '/dev/stdout', str(earth_radius), str(doppler_centroid)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             cwd=cwd, encoding='utf8')
        stdout_data, stderr_data = p.communicate(input=self.to_str())
        #print ('stdout_data', stdout_data)
        print (stderr_data)
        prm = PRM.from_str(stdout_data)
        if inplace:
            return self.set(prm)
        else:
            return prm

    def SAT_baseline(self, other, inplace=False):
        """
        SAT_baseline

        Usage: (two modes)

        mode 1:

            SAT_baseline PRM_master PRM_master

            This is used to compute height information
            (writes out height information for appending to PRM file)

        mode 2:

            SAT_baseline PRM_master PRM_aligned 

            PRM_master     PRM file for reference image
            PRM_aligned        PRM file of secondary image
            Please make sure the orbit file data is in PRM 
            Program runs through repeat orbit to find nearest point 
            to the start, center and end on the reference orbit
            (writes out parameters for appending to PRM file)
        """
        import os
        import subprocess

        if not isinstance(other, PRM):
            raise Exception('Argument "other" should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(self.to_str(), 'utf8'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(other.to_str(), 'utf8'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        argv = ['SAT_baseline', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}']
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=cwd, encoding='utf8')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        print (stderr_data)
        prm = PRM.from_str(stdout_data)
        if inplace:
            return self.set(prm)
        else:
            return prm

    """
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], precise=1)
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], precise=1, mode='-bos')
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], precise=1, mode='-bod')

    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], tofile='out.dat', precise=1)
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], tofile='outs.dat', precise=1, mode='-bos')
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], tofile='outd.dat', precise=1, mode='-bod')

    coords = prm.SAT_llt2rat(dem_data[:10], precise=1)
    [format(v, '.6f') for v in coords]
    """
    def SAT_llt2rat(self, coords=None, fromfile=None, tofile=None, precise=1, mode=None):
        """
         Usage: SAT_llt2rat master.PRM prec [-bo[s|d]] < inputfile > outputfile

             master.PRM   -  parameter file for master image and points to LED orbit file.
             precise      -  (0) standard back geocoding, (1) - polynomial refinenent (slower).
             inputfile    -  lon, lat, elevation [ASCII].
             outputfile   -  range, azimuth, elevation(ref to radius in PRM), lon, lat [ASCII default].
             -bos or -bod -  binary single or double precision output.
        """
        import numpy as np
        from io import StringIO, BytesIO
        import os
        import subprocess

        if mode is not None and not mode in ['-bos','-bod']:
            raise Exception('Supports only text (by default) and -bod/-bos for binary double and single precision output')

        if coords is not None and fromfile is None:
            #lon, lat, elevation = coords
            #data=f'{lon} {lat} {elevation}'
            buffer = BytesIO()
            # to produce the same result as original command returns
            # gmt grd2xyz --FORMAT_FLOAT_OUT=%lf topo/dem.grd -s
            np.savetxt(buffer, coords, delimiter=' ', fmt='%.6f')
            stdin_data = buffer.getvalue()
        elif coords is None and fromfile is not None:
            with open(fromfile, 'rb') as f:
                stdin_data = f.read()
        else:
            raise Exception('Should be defined data source as coordinates triplet (coords) or as file (fromfile)')

        # TBD: use np.array_split() and [x for x in a if x.size > 0] for chunking processing
        pipe = os.pipe()
        os.write(pipe[1], bytearray(self.to_str(), 'ascii'))
        os.close(pipe[1])
        #print ('descriptor', str(pipe[0]))

        argv = ['SAT_llt2rat', f'/dev/fd/{pipe[0]}', str(precise)]
        # set binary format mode
        if mode is not None:
            argv.append(mode)
        #print (argv)
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe[0]],
                             cwd=cwd, bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=stdin_data)

        stderr_data = stderr_data.decode('ascii')
        if stderr_data.startswith('interpolation point outside of data constraints'):
            print ('Error: SAT_llt2rat processing stopped due to invalid coordinates for one of input points')
            return None
        if stderr_data is not None and len(stderr_data):
            print (stderr_data)

        if tofile is not None:
            with open(tofile, 'wb') as f:
                f.write(stdout_data)
        else:
            if mode == '-bos':
                out = (np.frombuffer(stdout_data, dtype=np.dtype(np.float32)))
            if mode == '-bod':
                out = (np.frombuffer(stdout_data, dtype=np.dtype(np.float64)))
            else:
                out = np.fromstring(stdout_data, dtype=float, sep=' ')
            return out if out.size==5 else out.reshape(-1,5)

    def resamp(self, alignedPRM, alignedSLC_tofile, interp):
        """
        Usage: resamp master.PRM aligned.PRM new_aligned.PRM new_aligned.SLC interp
        master.PRM         - PRM for master imagea
        aligned.PRM        - PRM for aligned image
        new_aligned.PRM    - PRM for aligned aligned image
        new_aligned.SLC    - SLC for aligned aligned image
        interp             - interpolation method: 1-nearest; 2-bilinear; 3-biquadratic; 4-bisinc

        master  = PRM.from_file('master.PRM')
        aligned = PRM.from_file('aligned.PRM')
        master.resamp(aligned,'new_aligned.SLC',interp).to_file('new_aligned.PRM')
        """
        import os
        import subprocess

        if not isinstance(alignedPRM, PRM):
            raise Exception('Argument should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(alignedPRM.to_str(), 'utf8'))
        os.close(pipe1[1])
        #print ('descriptor 1', pipe1[0])

        pipe2 = os.pipe()
        #print ('descriptor 2', pipe2[1])

        # Usage: resamp master.PRM aligned.PRM new_aligned.PRM new_aligned.SLC interp
        #
        #cmd = f'resamp /dev/stdin /dev/fd/{pipe1[0]} /dev/fd/{pipe2[1]} /dev/stdout {intrp} | sponge ___'
        #cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        #p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        #                     stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[1]],
        #                     cwd=cwd, encoding='utf8', shell=True)
        argv = ['resamp', f'/dev/stdin', f'/dev/fd/{pipe1[0]}',
                f'/dev/fd/{pipe2[1]}', '/dev/stdout', str(interp)]
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[1]],
                             cwd=cwd)
        stdout_data, stderr_data = p.communicate(input=bytearray(self.to_str(), 'utf8'))

        # print errors and notifications
        print (stderr_data.decode('utf8'))

        # save big SLC binary file
        with open(alignedSLC_tofile, 'wb') as f:
            f.write(stdout_data)

        # PRM file should be small text
        data = os.read(pipe2[0],int(10e6)).decode('utf8')
        return PRM.from_str(data)

    # PRM.fitoffset(3, 3, offset_dat)
    # PRM.fitoffset(3, 3, matrix_fromfile='raw/offset.dat')
    @staticmethod
    def fitoffset(rank_rng, rank_azi, matrix=None, matrix_fromfile=None, SNR=20):
        import numpy as np
        import math

        """
        Usage: fitoffset.csh npar_rng npar_azi xcorr.dat [SNR]
 
            rank_rng    - number of parameters to fit in range 
            rank_azi    - number of parameters to fit in azimuth 
            matrix      - files of range and azimuth offset estimates 
            SNR         - optional SNR cutoff (default 20)
 
        Example: fitoffset.csh 3 3 freq_xcorr.dat 
        """
        if (matrix is None and matrix_fromfile is None) or (matrix is not None and matrix_fromfile is not None):
            raise Exception('One and only one argument matrix or matrix_fromfile should be defined')
        if matrix_fromfile is not None:
            matrix = np.genfromtxt(matrix_fromfile)

        #  first extract the range and azimuth data
        rng = matrix[np.where(matrix[:,4]>SNR)][:,[0,2,1]]
        azi = matrix[np.where(matrix[:,4]>SNR)][:,[0,2,3]]

        # make sure there are enough points remaining
        if rng.shape[0] < 8:
            raise Exception(f'FAILED - not enough points to estimate parameters, try lower SNR ({rng.shape[0]} < 8)')

        rng_coef = PRM.GMT_trend2d(rng, rank_rng)
        azi_coef = PRM.GMT_trend2d(azi, rank_azi)

        # range and azimuth data ranges
        scale_coef = [np.min(rng[:,0]), np.max(rng[:,0]), np.min(rng[:,1]), np.max(rng[:,1])]

        rng_coef += scale_coef
        azi_coef += scale_coef
        #print ('rng_coef, rng_coef)
        #print ('azi_coef, azi_coef)

        # now convert to range coefficients
        rshift = rng_coef[0] - rng_coef[1]*(rng_coef[4]+rng_coef[3])/(rng_coef[4]-rng_coef[3]) \
            - rng_coef[2]*(rng_coef[6]+rng_coef[5])/(rng_coef[6]-rng_coef[5])
        # now convert to azimuth coefficients
        ashift = azi_coef[0] - azi_coef[1]*(azi_coef[4]+azi_coef[3])/(azi_coef[4]-azi_coef[3]) \
            - azi_coef[2]*(azi_coef[6]+azi_coef[5])/(azi_coef[6]-azi_coef[5])
        #print ('rshift', rshift, 'ashift', ashift)

        # note: Python x % y expression and nympy results are different to C, use math function
        # use 'g' format for float values as in original GMTSAR codes to easy compare results
        prm = PRM().set(gformat=True,
                        rshift     =int(rshift) if rshift>=0 else int(rshift)-1,
                        sub_int_r  =math.fmod(rshift, 1)  if rshift>=0 else math.fmod(rshift, 1) + 1,
                        stretch_r  =rng_coef[1]*2/(rng_coef[4]-rng_coef[3]),
                        a_stretch_r=rng_coef[2]*2/(rng_coef[6]-rng_coef[5]),
                        ashift     =int(ashift) if ashift>=0 else int(ashift)-1,
                        sub_int_a  =math.fmod(ashift, 1)  if ashift>=0 else math.fmod(ashift, 1) + 1,
                        stretch_a  =azi_coef[1]*2/(azi_coef[4]-azi_coef[3]),
                        a_stretch_a=azi_coef[2]*2/(azi_coef[6]-azi_coef[5]),
                       )

        return prm

    # create replacement for trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
    def topo_ra(self, dem, trans_dat_tofile=None, topo_ra_tofile=None):
        import numpy as np
        import xarray as xr
        from scipy.interpolate import griddata

        XMAX, yvalid, num_patch = self.get('num_rng_bins', 'num_valid_az', 'num_patches')
        YMAX = yvalid * num_patch
        #print ('XMAX', XMAX, 'YMAX', YMAX)
        # as defined in dem2topo_ra.csh for Sentinel-1
        azi_dec = 2
        rng_dec = 2
        print (f'Range and azimuth decimation: {rng_dec}/{azi_dec}')

        # use center pixel GMT registration mode 
        rngs = np.arange(1,XMAX+1,2)
        azis = np.arange(1,YMAX+1,2)
        grid_r, grid_a = np.meshgrid(rngs, azis)
        #print ('grid_r', grid_r.shape)

        # create replacement for trans.dat - file generated by llt_grid2rat (r a topo lon lat)"
        lats, lons, z = xr.broadcast(dem.lat, dem.lon, dem)
        dem_data = np.column_stack([lons.values.ravel(), lats.values.ravel(), z.values.ravel()])
        #print ('dem_data', dem_data.shape)
        if trans_dat_tofile is not None:
            self.SAT_llt2rat(dem_data, tofile=trans_dat_tofile, precise=1, mode='-bod')
            coords = np.fromfile(trans_dat_tofile, dtype=np.float64).reshape([-1,5])
        else:
            coords = self.SAT_llt2rat(dem_data, precise=1, mode='-bod')

        grid = griddata((coords[:,0], coords[:,1]), coords[:,2], (grid_r, grid_a), method='linear')
        topo = xr.DataArray(np.flipud(grid), coords={'y': azis, 'x': rngs}, name='z')

        if topo_ra_tofile is not None:
            # save to NetCDF file
            compression = dict(zlib=True, complevel=3, chunksizes=[512,512])
            topo.to_netcdf(topo_ra_tofile, encoding={'z': compression})

        return topo

    def diff(self, other, gformat=True):
        """
        Compare to other dataframe and return difference
        """
        import pandas as pd
        import numpy as np

        if not isinstance(other, PRM):
            raise Exception('Argument should be PRM class instance')

        df1 = self.df.copy()
        df2 = other.df.copy()

        if gformat:
            fmt = lambda v: format(v, 'g') if type(v) in [float, np.float16, np.float32, np.float64, np.float128] else v
            df1['value'] = [fmt(value) for value in df1['value']]
            df2['value'] = [fmt(value) for value in df2['value']]

        return pd.concat([df1, df2]).drop_duplicates(keep=False)

    # TODO: add topo_ra argument processing
    # two binary files real.grd and imag.grd will be created
    # TBD: update phasediff tool to allow output files basename argument
    def phasediff(self, other, topo_ra_fromfile=None, topo_ra=None):
        """
        phasediff [GMTSAR] - Compute phase difference of two images

        Usage: phasediff ref.PRM rep.PRM [-topo topo_ra.grd] [-model modelphase.grd]
            (topo_ra and model in GMT grd format)
        """
        import os
        import subprocess

        if not isinstance(other, PRM):
            raise Exception('Argument "other" should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(self.to_str(), 'utf8'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(other.to_str(), 'utf8'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        argv = ['phasediff', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}']
        if topo_ra_fromfile is not None:
            argv.append('-topo')
            argv.append(topo_ra_fromfile)
        #print ('argv', argv)
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=cwd, encoding='utf8')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        print (stderr_data)
        return

def main():
    import sys

    if len(sys.argv) <= 2 :
        print (f"Usage: {sys.argv[0]} file.PRM parameter1 value1 parameter2 value2 ...")
        exit(0)

    prm_filename = sys.argv[1]
    #print ('prm_filename', prm_filename)
    prm = PRM.from_file(prm_filename)
    pairs = dict(zip(sys.argv[2::2], sys.argv[3::2]))
    prm = prm.set(**pairs)
    prm.to_file(prm_filename)

if __name__ == "__main__":
    # execute only if run as a script
    main()

"""
PRM.from_file('S1_20150403_ALL_F1.PRM').set(nrows=-999).update()

prm1 = PRM.from_file('S1_20150403_ALL_F1.PRM')
prm2 = PRM.from_file('S1_20150403_ALL_F1.PRM').sel('nrows') + 1000
(prm1).get('nrows'), (prm2).get('nrows'), (prm1 + prm2).get('nrows'), (prm1).get('nrows')

prm1.set(prm2).update()

(prm1 + prm2).to_file('S1_20150403_ALL_F1.PRM')

(prm1).get('nrows')

prm = PRM.from_file(...)
prm.set(prm.calc_dop_orb(0,0)).update()
"""
