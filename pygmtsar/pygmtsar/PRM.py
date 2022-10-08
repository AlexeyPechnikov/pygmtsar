#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to read, write and update PRM files and calculate Doppler orbit by calc_dop_orb command line tool
#import pytest

class PRM:

    # NetCDF options
    netcdf_chunksize = 512
    netcdf_compression = dict(zlib=True, complevel=3, chunksizes=(netcdf_chunksize,netcdf_chunksize))
    netcdf_engine = 'h5netcdf'

    @staticmethod
    def gmtsar_sharedir():
        """
        Call GMTSAR tool to obtain GMTSAR shared directory location
        """
        import os
        import subprocess

        argv = ['gmtsar_sharedir.csh']
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_data, stderr_data = p.communicate()

        stderr_data = stderr_data.decode('utf8')
        if stderr_data is not None and len(stderr_data):
            print ('gmtsar_sharedir', stderr_data)
        data = stdout_data.decode('utf8').strip()

        return data
    
    # replacement for GMTSAR gaussians
    # gauss5x5 = np.genfromtxt('/usr/local/GMTSAR/share/gmtsar/filters/gauss5x5',skip_header=True)
    # gaussian_kernel(5,1) ~= gauss5x5
    @staticmethod
    def gaussian_kernel(size=(5,5), std=(1,1)):
        """Make 2D Gaussian kernel matrix"""
        import numpy as np
        from scipy import signal
        matrix1 = signal.gaussian(size[0], std=std[0]).reshape(size[0], 1)
        matrix2 = signal.gaussian(size[1], std=std[1]).reshape(size[1], 1)
        matrix2d = np.outer(matrix1, matrix2)
        return matrix2d

    @staticmethod
    def nanconvolve2d(data, kernel, threshold=1/3.):
        """
        Convolution using generic kernel on a 2D array with NaNs
        """
        import numpy as np
        import scipy.signal
        # np.complex128 includes np.float64 real and imagine part
        vals = data.astype(np.complex128)
        #vals[np.isnan(data)] = 0
        #vals[~np.isnan(data)] += 1j
        nanmask = np.isnan(data)
        vals[nanmask] = 0
        vals[~nanmask] += 1j

        conv = scipy.signal.convolve2d(vals, 
                                      kernel.astype(np.complex128), 
                                      mode='same', boundary='symm'
                                     )
        # suppress incorrect division warning
        np.seterr(invalid='ignore')
        return np.where(conv.imag >= threshold*np.sum(kernel), np.divide(conv.real, conv.imag), np.nan)

    @staticmethod
    def nanconvolve2d_gaussian(data, sigmas, truncate):
        """
        Convolution using Gaussian kernel on a 2D array with NaNs
        """
        import numpy as np
        from scipy.ndimage import gaussian_filter

        # np.complex128 includes np.float64 real and imagine part
        vals = data.astype(np.complex128)
        nanmask = np.isnan(data)
        vals[nanmask] = 0
        vals[~nanmask] += 1j

        conv = gaussian_filter(vals, sigma=sigmas, truncate=truncate)
    
        # suppress incorrect division warning
        np.seterr(invalid='ignore')
        return np.divide(conv.real, conv.imag)

    @staticmethod
    def nearest_grid(in_grid, search_radius_pixels=300):
        """
        Pixel-based Nearest Neighbour interpolation
        """
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np
        ys, xs = np.meshgrid(range(in_grid.shape[1]), range(in_grid.shape[0]))
        ys = ys.reshape(-1)
        xs = xs.reshape(-1)
        if isinstance(in_grid, xr.DataArray):
            zs = in_grid.values.reshape(-1)
        else:
            zs = in_grid.reshape(-1)
        mask = np.where(~np.isnan(zs))
        # on regular source grid some options should be redefined for better performance
        tree = cKDTree(np.column_stack((ys[mask],xs[mask])), compact_nodes=False, balanced_tree=False)
        # use distance_limit
        d, inds = tree.query(np.column_stack((ys,xs)), k = 1, distance_upper_bound=search_radius_pixels, workers=8)
        # replace not available indexes by zero (see distance_upper_bound)
        fakeinds = np.where(~np.isinf(d), inds, 0)
        # produce the same output array as dataset to be able to add global attributes
        values = np.where(~np.isinf(d), zs[mask][fakeinds], np.nan).reshape(in_grid.shape)
        if isinstance(in_grid, xr.DataArray):
            return xr.DataArray(values, coords=in_grid.coords, name=in_grid.name)
        return values

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
                in [float, np.float16, np.float32, np.float64] else value
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
            #else:
            #    print (f'Safe mode: remain old PRM file {self.filename} and save new one {name}')
            # will be saved later
            self.filename = name

            input_file0 = os.path.join(dirname0, os.path.split(self.get('input_file'))[-1])
            input_ext = os.path.splitext(input_file0)[1]
            input_file = basename + input_ext
            #print (f'PRM change input_file attribute {input_file0} -> {input_file}')
            self.set(input_file = f'{shortname}{input_ext}')
            if os.path.isfile(input_file0) and not input_file == input_file0:
                if not safe:
                    #print (f'Rename input_file {input_file0} -> {input_file}')
                    if not debug:
                        os.rename(input_file0, input_file)
                else:
                    #print (f'Safe mode: copy input_file {input_file0} -> {input_file}')
                    if not debug:
                        shutil.copy2(input_file0, input_file, follow_symlinks=True)

            SLC_file0 = os.path.join(dirname0, os.path.split(self.get('SLC_file'))[-1])
            SLC_file = basename + '.SLC'
            #print (f'PRM change SLC_file attribute {SLC_file0} -> {SLC_file}')
            self.set(SLC_file = f'{shortname}.SLC')
            if os.path.isfile(SLC_file0) and not SLC_file == SLC_file0:
                #print (f'Rename SLC_file {SLC_file0} -> {SLC_file}')
                if not debug:
                    os.rename(SLC_file0, SLC_file)

            led_file0 = os.path.join(dirname0, os.path.split(self.get('led_file'))[-1])
            led_file = basename + '.LED'
            #print (f'PRM change LED_file attribute {led_file0} -> {led_file}')
            self.set(led_file = f'{shortname}.LED')
            if os.path.isfile(led_file0) and not led_file == led_file0:
                if not safe:
                    #print (f'Rename LED_file {led_file0} -> {led_file}')
                    if not debug:
                        os.rename(led_file0, led_file)
                else:
                    #print (f'Safe mode: copy LED_file {led_file0} -> {led_file}')
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

    def __sub__(self, other):
        import pandas as pd
        if isinstance(other, PRM):
            prm = pd.concat([self.df, other.df])
            # drop duplicates
            prm = prm.groupby(prm.index).last()
        else:
            prm = self.df - other
        return PRM(prm)

    def get(self, *args):
        out = [self.df.loc[[key]].iloc[0].values[0] for key in args]
        if len(out) == 1:
            return out[0]
        return out

    def shift_atime(self, lines, inplace=False):
        """
        Shift time in azimuth by a number of lines
        """
        prm = self.sel('clock_start','clock_stop','SC_clock_start','SC_clock_stop') + lines/self.get('PRF')/86400.0
        if inplace:
            return self.set(prm)
        else:
            return prm
    
    def calc_dop_orb(self, earth_radius=0, doppler_centroid=0, inplace=False, debug=False):
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
                             cwd=cwd, encoding='ascii')
        stdout_data, stderr_data = p.communicate(input=self.to_str())
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: calc_dop_orb', stderr_data)
        prm = PRM.from_str(stdout_data)
        if inplace:
            return self.set(prm)
        else:
            return prm

    def SAT_baseline(self, other, tail=None, debug=False):
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
        os.write(pipe1[1], bytearray(self.to_str(), 'ascii'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(other.to_str(), 'ascii'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        argv = ['SAT_baseline', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}']
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=cwd, encoding='ascii')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: SAT_baseline', stderr_data)
        prm = PRM.from_str(stdout_data)
        # replacement for SAT_baseline $1 $2 | tail -n9
        if tail is not None:
            prm.df = prm.df.tail(tail)
        return prm

    """
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], precise=1)
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], precise=1, binary=True)

    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], tofile='out.dat', precise=1)
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], tofile='outd.dat', precise=1, binary=True)

    coords = prm.SAT_llt2rat(dem_data[:10], precise=1)
    [format(v, '.6f') for v in coords]
    """
    def SAT_llt2rat(self, coords=None, fromfile=None, tofile=None, precise=1, binary=False, debug=False):
        """
         Usage: SAT_llt2rat master.PRM prec [-bo[s|d]] < inputfile > outputfile

             master.PRM   -  parameter file for master image and points to LED orbit file.
             precise      -  (0) standard back geocoding, (1) - polynomial refinenent (slower).
             inputfile    -  lon, lat, elevation [ASCII].
             outputfile   -  range, azimuth, elevation(ref to radius in PRM), lon, lat [ASCII default].
             -bos or -bod -  binary single or double precision output.

             Note: -bos mode support deleted as obsolete
        """
        import numpy as np
        from io import StringIO, BytesIO
        import os
        import subprocess

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
        if binary:
            argv.append('-bod')
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
        if stderr_data is not None and len(stderr_data) and debug:
            print ('DEBUG: SAT_llt2rat', stderr_data)

        if tofile is not None:
            with open(tofile, 'wb') as f:
                f.write(stdout_data)
        else:
            if binary:
                out = (np.frombuffer(stdout_data, dtype=np.dtype(np.float64)))
            else:
                out = np.fromstring(stdout_data, dtype=float, sep=' ')
            return out if out.size==5 else out.reshape(-1,5)

    def resamp(self, alignedPRM, alignedSLC_tofile, interp, debug=False):
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
        os.write(pipe1[1], bytearray(alignedPRM.to_str(), 'ascii'))
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
        #                     cwd=cwd, encoding='ascii', shell=True)
        argv = ['resamp', f'/dev/stdin', f'/dev/fd/{pipe1[0]}',
                f'/dev/fd/{pipe2[1]}', '/dev/stdout', str(interp)]
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[1]],
                             cwd=cwd)
        stdout_data, stderr_data = p.communicate(input=bytearray(self.to_str(), 'ascii'))

        # print errors and notifications
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: resamp', stderr_data.decode('ascii'))

        # save big SLC binary file
        with open(alignedSLC_tofile, 'wb') as f:
            f.write(stdout_data)

        # PRM file should be small text
        data = os.read(pipe2[0],int(10e6)).decode('ascii')
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
            fmt = lambda v: format(v, 'g') if type(v) in [float, np.float16, np.float32, np.float64] else v
            df1['value'] = [fmt(value) for value in df1['value']]
            df2['value'] = [fmt(value) for value in df2['value']]

        return pd.concat([df1, df2]).drop_duplicates(keep=False)

    # TODO: add topo_ra argument processing
    # two binary files real.grd and imag.grd will be created
    # TBD: update phasediff tool to allow output files basename argument
    def phasediff(self, other, topo_ra_fromfile, imag_tofile, real_tofile, debug=False):
        """
        phasediff [GMTSAR] - Compute phase difference of two images

        Usage: phasediff ref.PRM rep.PRM [-topo topo_ra.grd] [-model modelphase.grd] [-imag imag.grd] [-real real.grd]
            (topo_ra and model in GMT grd format)
        """
        import os
        import subprocess

        if not isinstance(other, PRM):
            raise Exception('Argument "other" should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(self.to_str(), 'ascii'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(other.to_str(), 'ascii'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        argv = ['phasediff', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}',
                '-imag', os.path.relpath(imag_tofile, cwd),
                '-real', os.path.relpath(real_tofile, cwd)]
        if topo_ra_fromfile is not None:
            argv.append('-topo')
            argv.append(os.path.relpath(topo_ra_fromfile,cwd))
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=cwd, encoding='ascii')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: phasediff', stderr_data)
        return

    # Usage: make_gaussian_filter name_of_PRM_file RNG_DEC AZI_DEC WAVELENGTH(m)
    # make_gaussian_filter S1_20150121_ALL_F1.PRM 2 1 400
    def make_gaussian_filter(self, range_dec, azi_dec, wavelength=200, debug=False):
        """
        make_gaussian_filter

        Usage: make_gaussian_filter name_of_PRM_file RNG_DEC AZI_DEC WAVELENGTH(m) [output_filename]
        Example: make_gaussian_filter IMG-HH-ALPSRP211830620-H1.0__A.PRM 2 4 200
        Output: gauss_WAVELENGTH or output_filename
        """
        import os
        import subprocess
        import numpy as np

        pipe1 = os.pipe()
        #print ('descriptor 1', pipe1[1])

        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        argv = ['make_gaussian_filter', '/dev/stdin',
                str(range_dec), str(azi_dec), str(wavelength), f'/dev/fd/{pipe1[1]}']
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[1]],
                             cwd=cwd, encoding='ascii')
        stdout_data, stderr_data = p.communicate(self.to_str())
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: make_gaussian_filter', stderr_data)
            print ('DEBUG: make_gaussian_filter', stdout_data)

        data = os.read(pipe1[0],int(10e6)).decode('ascii')
        return [np.fromstring(stdout_data, dtype=int, sep=' '), data]

    # the command line tool supports only GMT binary format input (=bf) while NetCDF output is possible too
    def conv(self, idec, jdec, output_file, filter_file=None, filter_string=None, input_file=None, debug=False):
        """
        conv [GMTSAR] - 2-D image convolution

        Usage: conv idec jdec filter_file input output
           idec           - row decimation factor
           jdec           - column decimation factor
           filter_file    - eg. filters/gauss17x5
           input          - name of file to be filtered (I*2 or R*4)
           output         - name of filtered output file (R*4 only)

           examples:
           conv 4 2 filters/gauss9x5 IMG-HH-ALPSRP109430660-H1.0__A.PRM test.grd
           (makes and filters amplitude file from an SLC-file)

           conv 4 2 filters/gauss5x5 infile.grd outfile.grd
           (filters a float file)
        """
        import os
        import subprocess

        if filter_file is None and filter_string is None:
            raise Exception('Should be defined argument filter_file or filter_string')

        if input_file is None and self.filename is None:
            raise Exception('PRM object should be created from file or defined argument input_file')

        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        if input_file is None:
            input_filename = os.path.basename(self.filename)
        else:
            input_filename = os.path.relpath(input_file, cwd)

        # open descriptors
        fds = []
        if filter_string is not None:
            pipe2 = os.pipe()
            os.write(pipe2[1], bytearray(filter_string, 'ascii'))
            os.close(pipe2[1])
            input_filter_filename = f'/dev/fd/{pipe2[0]}'
            fds.append(pipe2[0])
            #print ('descriptor 2', str(pipe2[0]))
        else:
            input_filter_filename = os.path.relpath(filter_file, cwd)

        if output_file is None:
            output_filename = '/dev/stdout'
        else:
            output_filename = os.path.relpath(output_file, cwd)

        argv = ['conv', str(idec), str(jdec), input_filter_filename, input_filename, output_filename]
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=fds,
                             cwd=cwd)
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', len(stdout_data))
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: conv', stderr_data.decode('ascii'))
        return

    # actually, -alpha command line argument is useless
    # if amp files are defined set alpha = -1 to denote coherence-based alpha
    def phasefilt(self, imag_fromfile, real_fromfile, amp1_fromfile, amp2_fromfile,
                  phasefilt_tofile, corrfilt_tofile, psize=32, debug=False):
        """
        phasefilt [GMTSAR] - Apply adaptive non-linear phase filter

        USAGE:
            phasefilt -imag imag.grd -real real.grd [-alpha alpha][-psize size][-amp1 amp1.grd -amp2 amp2.grd][-diff][-v]
             applies Goldstein adaptive filter to phase [output: filtphase.grd]
             or applies modified Goldstein adaptive filter to phase [output: filtphase.grd, corrfilt.grd]
            -imag [required] GMT format file of imaginary component
            -real [required] GMT format file of real component
            -alpha  exponent for filter - usually between 0.0 and 1.5 (0.0 should not filter).
                    default: 0.5 [Goldstein filter] (anything above 1.0 may be excessive)
                    alpha < 0 will set alpha = (1 - coherence) [modified Goldstein]
            -psize patch size for filtering. Must be power of two.
                    default: 32
            -amp1 GMT format file of amplitude image of image 1. Needed (and applies) modified filter.
            -amp2 GMT format file of amplitude image of image 2. Needed (and applies) modified filter.
        """
        import os
        import subprocess

        if corrfilt_tofile is None and alpha<0:
            raise Exception('For argument alpha<0 shoukd be defined argument corrfilt_tofile')

        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        # -v - only args after verbose (if set) will be echoed
        argv = ['phasefilt',
                '-imag', os.path.relpath(imag_fromfile, cwd),
                '-real', os.path.relpath(real_fromfile, cwd),
                '-amp1', os.path.relpath(amp1_fromfile, cwd),
                '-amp2', os.path.relpath(amp2_fromfile, cwd),
                '-phasefilt', os.path.relpath(phasefilt_tofile, cwd),
                '-corrfilt', os.path.relpath(corrfilt_tofile, cwd),
                '-psize', str(psize)]
        #print ('argv', argv)

        p = subprocess.Popen(argv, stderr=subprocess.PIPE, cwd=cwd)
        stderr_data = p.communicate()[1]
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: phasefilt', stderr_data.decode('ascii'))
        return

    # see about correlation filter
    # https://github.com/gmtsar/gmtsar/issues/86
    # use_boxcar_filter=True for ISCE-type boxcar and multilook filter
    def intf(self, other, basedir, topo_ra_fromfile, basename=None, wavelength=200, psize=32,
             use_boxcar_filter=False, func=None, chunks='auto'):
        import os
        import numpy as np
        import xarray as xr
        #from scipy import signal
        import dask_image.ndfilters
        import dask.array

        # constant from GMTSAR code
        thresh = 5.e-21

        if not isinstance(other, PRM):
            raise Exception('Argument "other" should be PRM class instance')

        # define basename from two PRM names
        if basename is None:
            # S1_20171111_ALL_F3.PRM -> F3
            subswath = os.path.basename(self.filename)[16:18]
            date1    = os.path.basename(self.filename)[3:11]
            date2    = os.path.basename(other.filename)[3:11]
            basename = f'{subswath}_{date1}_{date2}_'
            #print ('basename', basename)

        # make full file name
        fullname = lambda name: os.path.join(basedir, basename + name)
        gmtsar_sharedir = self.gmtsar_sharedir()

        filename_fill_3x3  = os.path.join(gmtsar_sharedir, 'filters', 'fill.3x3')
        filename_boxcar3x5 = os.path.join(gmtsar_sharedir, 'filters', 'boxcar.3x5')
        filename_gauss5x5  = os.path.join(gmtsar_sharedir, 'filters', 'gauss5x5')
        
        #!conv 1 1 /usr/local/GMTSAR/share/gmtsar/filters/fill.3x3 raw/tmp2.nc raw/corr.nc
        fill_3x3 = np.genfromtxt(filename_fill_3x3, skip_header=1)
        gauss_dec, gauss_string = self.make_gaussian_filter(2, 1, wavelength=wavelength)
        #print (gauss_matrix_astext)

        # prepare PRMs for the calculation below
        other.set(self.SAT_baseline(other, tail=9))
        self.set(self.SAT_baseline(self).sel('SC_height','SC_height_start','SC_height_end'))

        # for topo_ra use relative path from PRM files directory
        # use imag.grd=bf for GMT native, C-binary format
        # os.path.join(basedir, f'{subswath}_topo_ra.grd')
        self.phasediff(other, topo_ra_fromfile=topo_ra_fromfile,
                       imag_tofile=fullname('imag.grd=bf'),
                       real_tofile=fullname('real.grd=bf'))

        # filtering interferogram
        if not use_boxcar_filter:
            # use default GMTSAR filter
            # making amplitudes
            self.conv(1, 2, filter_file = filename_gauss5x5,
                      output_file=fullname('amp1_tmp.grd=bf'))
            self.conv(gauss_dec[0], gauss_dec[1], filter_string=gauss_string,
                      input_file=fullname('amp1_tmp.grd=bf'),
                      output_file=fullname('amp1.grd'))

            other.conv(1, 2, filter_file = filename_gauss5x5,
                       output_file=fullname('amp2_tmp.grd=bf'))
            other.conv(gauss_dec[0], gauss_dec[1], filter_string=gauss_string,
                       input_file=fullname('amp2_tmp.grd=bf'),
                       output_file=fullname('amp2.grd'))

            # filtering interferogram
            self.conv(1, 2, filter_file=filename_gauss5x5,
                      input_file=fullname('real.grd=bf'),
                      output_file=fullname('real_tmp.grd=bf'))
            other.conv(gauss_dec[0], gauss_dec[1], filter_string=gauss_string,
                       input_file=fullname('real_tmp.grd=bf'),
                       output_file=fullname('realfilt.grd'))
            self.conv(1, 2, filter_file=filename_gauss5x5,
                      input_file=fullname('imag.grd=bf'),
                      output_file=fullname('imag_tmp.grd=bf'))
            other.conv(gauss_dec[0], gauss_dec[1], filter_string=gauss_string,
                       input_file=fullname('imag_tmp.grd=bf'),
                       output_file=fullname('imagfilt.grd'))
        else:
            # use ISCE-type boxcar and multilook filter
            # 3 range and 5 azimuth looks
            # making amplitudes
            self.conv(5, 3, filter_file = filename_boxcar3x5,
                      output_file=fullname('amp1.grd'))
            other.conv(5, 3, filter_file = filename_boxcar3x5,
                       output_file=fullname('amp2.grd'))

            self.conv(5, 3, filter_file=filename_boxcar3x5,
                      input_file=fullname('real.grd=bf'),
                      output_file=fullname('realfilt.grd'))
            self.conv(5, 3, filter_file=filename_boxcar3x5,
                      input_file=fullname('imag.grd=bf'),
                      output_file=fullname('imagfilt.grd'))

        # filtering phase
        self.phasefilt(imag_fromfile=fullname('imagfilt.grd'),
                       real_fromfile=fullname('realfilt.grd'),
                       amp1_fromfile=fullname('amp1.grd'),
                       amp2_fromfile=fullname('amp2.grd'),
                       phasefilt_tofile=fullname('phasefilt_phase.grd'),
                       corrfilt_tofile=fullname('phasefilt_corr.grd'),
                       psize=psize)

        # Python post-processing
        # we need to flip vertically results from the command line tools
        realfilt = xr.open_dataarray(fullname('realfilt.grd'), engine=self.netcdf_engine, chunks=chunks)
        imagfilt = xr.open_dataarray(fullname('imagfilt.grd'), engine=self.netcdf_engine, chunks=chunks)
        amp = np.sqrt(realfilt**2 + imagfilt**2)

        amp1 = xr.open_dataarray(fullname('amp1.grd'), engine=self.netcdf_engine, chunks=chunks)
        amp2 = xr.open_dataarray(fullname('amp2.grd'), engine=self.netcdf_engine, chunks=chunks)

        # use the same coordinates for all output grids
        # use .values to remove existing attributes from the axes
        coords = {'y': amp.y.values, 'x': amp.x.values}

        # making correlation
        tmp = amp1 * amp2
        mask = xr.where(tmp >= thresh, 1, np.nan)
        tmp2 = mask * (amp/np.sqrt(tmp))

        #conv = signal.convolve2d(tmp2, fill_3x3/fill_3x3.sum(), mode='same', boundary='symm')
        # use dask rolling window for the same convolution - 1 border pixel is NaN here
        #kernel = xr.DataArray(fill_3x3, dims=['i', 'j'])/fill_3x3.sum()
        #conv = tmp2.rolling(y=3, x=3, center={'y': True, 'x': True}).construct(lat='j', lon='i').dot(kernel)
        # use dask_image package
        kernel = xr.DataArray(fill_3x3, dims=['y', 'x'])/fill_3x3.sum()
        conv = dask_image.ndfilters.convolve(tmp2.data, kernel.data, mode='reflect')
        
        # wrap dask or numpy array to dataarray
        corr = xr.DataArray(dask.array.flipud(conv.astype(np.float32)), coords, name='z')
        if func is not None:
            corr = func(corr)
        if os.path.exists(fullname('corr.grd')):
            os.remove(fullname('corr.grd'))
        corr.to_netcdf(fullname('corr.grd'), encoding={'z': self.netcdf_compression}, engine=self.netcdf_engine)

        # make the Werner/Goldstein filtered phase
        phasefilt_phase = xr.open_dataarray(fullname('phasefilt_phase.grd'), engine=self.netcdf_engine, chunks=chunks)
        phasefilt = phasefilt_phase * mask
        phasefilt = xr.DataArray(dask.array.flipud(phasefilt.astype(np.float32)), coords, name='z')
        if func is not None:
            phasefilt = func(phasefilt)
        if os.path.exists(fullname('phasefilt.grd')):
            os.remove(fullname('phasefilt.grd'))
        phasefilt.to_netcdf(fullname('phasefilt.grd'), encoding={'z': self.netcdf_compression}, engine=self.netcdf_engine)

        # cleanup
        for name in ['amp1_tmp.grd', 'amp2_tmp.grd', 'amp1.grd', 'amp2.grd',
                     'real.grd', 'real_tmp.grd', 'realfilt.grd',
                     'imag.grd', 'imag_tmp.grd', 'imagfilt.grd',
                     'phasefilt_phase.grd', 'phasefilt_corr.grd']:
            filename = fullname(name)
            if not os.path.exists(filename):
                continue
            os.remove(filename)

        return

    def SAT_look(self, coords=None, fromfile=None, tofile=None, binary=False, debug=False):
        """
        Usage: SAT_look master.PRM [-bo[s|d]] < inputfile > outputfile

            master.PRM   -  parameter file for master image and points to LED orbit file
            inputfile    -  lon, lat, elevation [ASCII]
            outputfile   -  lon, lat, elevation look_E look_N look_U [ASCII default]
            -bos or -bod -  binary single or double precision output

                Note: -bos mode does not work

        example: SAT_look master.PRM < topo.llt > topo.lltn

        Note that the output elevation is the one above reference radius specified in the PRM file
        """
        import numpy as np
        from io import StringIO, BytesIO
        import os
        import subprocess

        if coords is not None and fromfile is None:
            #lon, lat, elevation = coords
            #data=f'{lon} {lat} {elevation}'
            buffer = BytesIO()
            np.savetxt(buffer, coords, delimiter=' ', fmt='%.6f')
            stdin_data = buffer.getvalue()
        elif coords is None and fromfile is not None:
            with open(fromfile, 'rb') as f:
                stdin_data = f.read()
        else:
            raise Exception('Should be defined data source as coordinates triplet (coords) or as file (fromfile)')

        pipe = os.pipe()
        os.write(pipe[1], bytearray(self.to_str(), 'ascii'))
        os.close(pipe[1])
        #print ('descriptor', str(pipe[0]))

        argv = ['SAT_look', f'/dev/fd/{pipe[0]}']
        # set binary format mode
        if binary:
            argv.append('-bod')
        #print (argv)
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe[0]],
                             cwd=cwd, bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=stdin_data)

        stderr_data = stderr_data.decode('ascii')
        if stderr_data is not None and len(stderr_data) and debug:
            print ('DEBUG: SAT_look', stderr_data)
            return None

        if tofile is not None:
            with open(tofile, 'wb') as f:
                f.write(stdout_data)
        else:
            if binary:
                out = (np.frombuffer(stdout_data, dtype=np.dtype(np.float64)))
            else:
                out = np.fromstring(stdout_data, dtype=float, sep=' ')
            return out if out.size==6 else out.reshape(-1,6)

    def pixel_size(self):
        """
        Calculate azimuth and range pixel size in meters
        Note: see make_gaussian_filter.c for the original code
        """
        import numpy as np
        from scipy.constants import speed_of_light

        # compute the range and azimuth pixel size
        RE, vel, ht, prf, fs, near_range, num_rng_bins = \
            self.get('earth_radius','SC_vel', 'SC_height', 'PRF', 'rng_samp_rate', 'near_range', 'num_rng_bins')
        #ER, vel, ht, prf, fs, near_range, num_rng_bins
        # real_vel/prf
        azi_px_size = vel / np.sqrt(1 + ht / RE) / prf
        rng_px_size = speed_of_light / fs / 2.0
        # compute the cosine of the looking angle and the surface deviate angle
        a = ht + RE
        far_range = near_range + rng_px_size * num_rng_bins
        rng = (near_range + far_range) / 2.0
        cost = (np.power(a, 2.0) + np.power(rng, 2.0) - np.power(RE, 2.0)) / 2.0 / a / rng
        cosa = (np.power(a, 2.0) + np.power(RE, 2.0) - np.power(rng, 2.0)) / 2.0 / a / RE
        # compute the ground range pixel size
        rng_px_size = rng_px_size / np.sin(np.arccos(cost) + np.arccos(cosa))
        # ground spacing in meters
        return (azi_px_size, rng_px_size)

    # TODO: use PRM parameters to define config parameters
    def snaphu_config(self, defomax=0, **kwargs):
        import os
        # we already use joblib everywhere
        import joblib

        tiledir = os.path.splitext(self.filename)[0]
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
    DEFOMAX_CYCLE  {defomax}
    TILEDIR        {tiledir}_snaphu_tiledir
    NPROC          {n_jobs}
    """
        conf_custom = '# custom config\n'
        for key, value in kwargs.items():
            if isinstance(value, bool):
                value = 'TRUE' if value else 'FALSE'
            conf_custom += f'    {key} {value}\n'
        return conf_basic + conf_custom

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
