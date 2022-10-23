#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .datagrid import datagrid
from .PRM_gmtsar import PRM_gmtsar

class PRM(datagrid, PRM_gmtsar):

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

    # see about correlation filter
    # https://github.com/gmtsar/gmtsar/issues/86
    # use_boxcar_filter=True for ISCE-type boxcar and multilook filter
    def intf(self, other, basedir, topo_ra_fromfile, basename=None, wavelength=200, psize=32,
             use_boxcar_filter=False, func=None, chunks='auto', debug=False):
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
                       real_tofile=fullname('real.grd=bf'),
                       debug=debug)

        # filtering interferogram
        if not use_boxcar_filter:
            # use default GMTSAR filter
            # making amplitudes
            self.conv(1, 2, filter_file = filename_gauss5x5,
                      output_file=fullname('amp1_tmp.grd=bf'),
                      debug=debug)
            self.conv(gauss_dec[0], gauss_dec[1], filter_string=gauss_string,
                      input_file=fullname('amp1_tmp.grd=bf'),
                      output_file=fullname('amp1.grd'),
                      debug=debug)

            other.conv(1, 2, filter_file = filename_gauss5x5,
                       output_file=fullname('amp2_tmp.grd=bf'),
                       debug=debug)
            other.conv(gauss_dec[0], gauss_dec[1], filter_string=gauss_string,
                       input_file=fullname('amp2_tmp.grd=bf'),
                       output_file=fullname('amp2.grd'),
                       debug=debug)

            # filtering interferogram
            self.conv(1, 2, filter_file=filename_gauss5x5,
                      input_file=fullname('real.grd=bf'),
                      output_file=fullname('real_tmp.grd=bf'),
                      debug=debug)
            other.conv(gauss_dec[0], gauss_dec[1], filter_string=gauss_string,
                       input_file=fullname('real_tmp.grd=bf'),
                       output_file=fullname('realfilt.grd'),
                       debug=debug)
            self.conv(1, 2, filter_file=filename_gauss5x5,
                      input_file=fullname('imag.grd=bf'),
                      output_file=fullname('imag_tmp.grd=bf'),
                      debug=debug)
            other.conv(gauss_dec[0], gauss_dec[1], filter_string=gauss_string,
                       input_file=fullname('imag_tmp.grd=bf'),
                       output_file=fullname('imagfilt.grd'),
                       debug=debug)
        else:
            # use ISCE-type boxcar and multilook filter
            # 3 range and 5 azimuth looks
            # making amplitudes
            self.conv(5, 3, filter_file = filename_boxcar3x5,
                      output_file=fullname('amp1.grd'),
                      debug=debug)
            other.conv(5, 3, filter_file = filename_boxcar3x5,
                       output_file=fullname('amp2.grd'),
                       debug=debug)

            self.conv(5, 3, filter_file=filename_boxcar3x5,
                      input_file=fullname('real.grd=bf'),
                      output_file=fullname('realfilt.grd'),
                      debug=debug)
            self.conv(5, 3, filter_file=filename_boxcar3x5,
                      input_file=fullname('imag.grd=bf'),
                      output_file=fullname('imagfilt.grd'),
                      debug=debug)

        # filtering phase
        self.phasefilt(imag_fromfile=fullname('imagfilt.grd'),
                       real_fromfile=fullname('realfilt.grd'),
                       amp1_fromfile=fullname('amp1.grd'),
                       amp2_fromfile=fullname('amp2.grd'),
                       phasefilt_tofile=fullname('phasefilt_phase.grd'),
                       corrfilt_tofile=fullname('phasefilt_corr.grd'),
                       psize=psize,
                       debug=debug)

        # Python post-processing
        # we need to flip vertically results from the command line tools
        realfilt = xr.open_dataarray(fullname('realfilt.grd'), engine=self.engine, chunks=chunks)
        imagfilt = xr.open_dataarray(fullname('imagfilt.grd'), engine=self.engine, chunks=chunks)
        amp = np.sqrt(realfilt**2 + imagfilt**2)

        amp1 = xr.open_dataarray(fullname('amp1.grd'), engine=self.engine, chunks=chunks)
        amp2 = xr.open_dataarray(fullname('amp2.grd'), engine=self.engine, chunks=chunks)

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
        corr.to_netcdf(fullname('corr.grd'), encoding={'z': self.compression}, engine=self.engine)

        # make the Werner/Goldstein filtered phase
        phasefilt_phase = xr.open_dataarray(fullname('phasefilt_phase.grd'), engine=self.engine, chunks=chunks)
        phasefilt = phasefilt_phase * mask
        phasefilt = xr.DataArray(dask.array.flipud(phasefilt.astype(np.float32)), coords, name='z')
        if func is not None:
            phasefilt = func(phasefilt)
        if os.path.exists(fullname('phasefilt.grd')):
            os.remove(fullname('phasefilt.grd'))
        phasefilt.to_netcdf(fullname('phasefilt.grd'), encoding={'z': self.compression}, engine=self.engine)

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
