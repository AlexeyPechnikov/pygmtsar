# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .datagrid import datagrid
from .PRM_gmtsar import PRM_gmtsar

class PRM(datagrid, PRM_gmtsar):

    # my replacement function for GMT based robust 2D trend coefficient calculations:
    # gmt trend2d r.xyz -Fxyzmw -N1r -V
    # gmt trend2d r.xyz -Fxyzmw -N2r -V
    # gmt trend2d r.xyz -Fxyzmw -N3r -V
    # https://github.com/GenericMappingTools/gmt/blob/master/src/trend2d.c#L719-L744
    # 3 model parameters
    # rank = 3 => nu = size-3
    @staticmethod
    def robust_trend2d(data, rank):
        """
        Perform robust linear regression to estimate the trend in 2D data.

        Parameters
        ----------
        data : numpy.ndarray
            Array containing the input data. The shape of the array should be (N, 3), where N is the number of data points.
            The first column represents the x-coordinates, the second column represents the y-coordinates (if rank is 3),
            and the third column represents the z-values.
        rank : int
            Number of model parameters to fit. Should be 1, 2, or 3. If rank is 1, the function fits a constant trend.
            If rank is 2, it fits a linear trend. If rank is 3, it fits a planar trend.

        Returns
        -------
        numpy.ndarray
            Array containing the estimated trend coefficients. The length of the array depends on the specified rank.
            For rank 1, the array will contain a single value (intercept).
            For rank 2, the array will contain two values (intercept and slope).
            For rank 3, the array will contain three values (intercept, slope_x, slope_y).

        Raises
        ------
        Exception
            If the specified rank is not 1, 2, or 3.

        Notes
        -----
        The function performs robust linear regression using the M-estimator technique. It iteratively fits a linear model
        and updates the weights based on the residuals until convergence. The weights are adjusted using Tukey's bisquare
        weights to downweight outliers.

        References
        ----------
        - Rousseeuw, P. J. (1984). Least median of squares regression. Journal of the American statistical Association, 79(388), 871-880.

        - Huber, P. J. (1973). Robust regression: asymptotics, conjectures and Monte Carlo. The Annals of Statistics, 1(5), 799-821.
        """
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

#     # standalone function is compatible but it is too slow while it should not be a problem
#     @staticmethod
#     def robust_trend2d(data, rank):
#         import numpy as np
#         import statsmodels.api as sm
# 
#         if rank not in [1, 2, 3]:
#             raise Exception('Number of model parameters "rank" should be 1, 2, or 3')
# 
#         if rank in [2, 3]:
#             x = data[:, 0]
#             x = np.interp(x, (x.min(), x.max()), (-1, +1))
#         if rank == 3:
#             y = data[:, 1]
#             y = np.interp(y, (y.min(), y.max()), (-1, +1))
#         z = data[:, 2]
# 
#         if rank == 1:
#             X = sm.add_constant(np.ones(z.shape))
#         elif rank == 2:
#             X = sm.add_constant(x)
#         elif rank == 3:
#             X = sm.add_constant(np.column_stack((x, y)))
# 
#         # Hampel weighting function looks the most accurate for noisy data
#         rlm_model = sm.RLM(z, X, M=sm.robust.norms.Hampel())
#         rlm_results = rlm_model.fit()
# 
#         return rlm_results.params[:rank]

#     # calculate MSE (optional)
#     @staticmethod
#     def robust_trend2d_mse(data, coeffs, rank):
#         import numpy as np
# 
#         x = data[:, 0]
#         y = data[:, 1]
#         z_actual = data[:, 2]
# 
#         if rank == 1:
#             z_predicted = coeffs[0]
#         elif rank == 2:
#             z_predicted = coeffs[0] + coeffs[1] * x
#         elif rank == 3:
#             z_predicted = coeffs[0] + coeffs[1] * x + coeffs[2] * y
# 
#         mse = np.mean((z_actual - z_predicted) ** 2)
#         return mse

    @staticmethod
    def from_list(prm_list):
        """
        Convert a list of parameter and value pairs to a PRM object.

        Parameters
        ----------
        prm_list : list
            A list of PRM strings.

        Returns
        -------
        PRM
            A PRM object.
        """
        from io import StringIO
        prm = StringIO('\n'.join(prm_list))
        return PRM._from_io(prm)

    @staticmethod
    def from_str(prm_string):
        """
        Convert a string of parameter and value pairs to a PRM object.

        Parameters
        ----------
        prm_string : str
            A PRM string.

        Returns
        -------
        PRM
            A PRM object.
        """
        from io import StringIO
        prm = StringIO(prm_string)
        return PRM._from_io(prm)

    @staticmethod
    def from_file(prm_filename):
        """
        Convert a PRM file of parameter and value pairs to a PRM object.

        Parameters
        ----------
        prm_filename : str
            The filename of the PRM file.

        Returns
        -------
        PRM
            A PRM object.
        """
        #data = json.loads(document)
        prm = PRM._from_io(prm_filename)
        prm.filename = prm_filename
        return prm

    @staticmethod
    def _from_io(prm):
        """
        Read parameter and value pairs from IO stream to a PRM object.

        Parameters
        ----------
        prm : IO stream
            The IO stream.

        Returns
        -------
        PRM
            A PRM object.
        """
        import pandas as pd
        return PRM(pd.read_csv(prm, sep='\s+=\s+', header=None, names=['name', 'value'], engine='python').set_index('name')\
                    .applymap(lambda val : pd.to_numeric(val,errors='ignore')))

    def __init__(self, prm=None):
        """
        Initialize a PRM object.

        Parameters
        ----------
        prm : PRM or pd.DataFrame, optional
            The PRM object or DataFrame to initialize from. Default is None.

        Returns
        -------
        None
        """
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
        """
        Compare two PRM objects for equality.

        Parameters
        ----------
        other : PRM
            The other PRM object to compare with.

        Returns
        -------
        bool
            True if the PRM objects are equal, False otherwise.
        """
        return isinstance(self, PRM) and self.df == other.df

    def __str__(self):
        """
        Return a string representation of the PRM object.

        Returns
        -------
        str
            The string representation of the PRM object.
        """
        return self.to_str()

    def __repr__(self):
        """
        Return a string representation of the PRM object for debugging.

        Returns
        -------
        str
            The string representation of the PRM object. If the PRM object was created from a file, 
            the filename and the number of items in the DataFrame representation of the PRM object 
            are included in the string.
        """
        if self.filename:
            return 'Object %s (%s) %d items\n%r' % (self.__class__.__name__, self.filename, len(self.df), self.df)
        else:
            return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    # use 'g' format for Python and numpy float values
    def set(self, prm=None, gformat=False, **kwargs):
        """
        Set PRM values.

        Parameters
        ----------
        prm : PRM, optional
            The PRM object to set values from. Default is None.
        gformat : bool, optional
            Whether to use 'g' format for float values. Default is False.
        **kwargs
            Additional keyword arguments for setting individual values.

        Returns
        -------
        PRM
            The updated PRM object.
        """
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
        """
        Convert the PRM object to a DataFrame.

        Returns
        -------
        pd.DataFrame
            The DataFrame representation of the PRM object.
        """
        return self.df

    def to_file(self, prm):
        """
        Save the PRM object to a PRM file.

        Parameters
        ----------
        prm : str
            The filename of the PRM file to save to.

        Returns
        -------
        PRM
            The PRM object.
        """
        self._to_io(prm)
        return self

    #def update(self):
    #    if self.filename is None:
    #        raise Exception('PRM is not created from file, use to_file() method instead')
    #    return self._to_io(self.filename)
    def update(self, name=None, safe=False, debug=False):
        """
        Save PRM file to disk and rename all 3 linked files together: input, LED, SLC if "name" defined.
        If safe=True, save old (small) PRM and LED files and move only (large) SLC file.

        Parameters
        ----------
        name : str, optional
            The new name for the PRM file. Default is None.
        safe : bool, optional
            Whether to use safe mode. Default is False.
        debug : bool, optional
            Whether to enable debug mode. Default is False.

        Returns
        -------
        PRM
            The updated PRM object.
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
        """
        Convert the PRM object to a string.

        Returns
        -------
        str
            The PRM string.
        """
        return self._to_io()

    def _to_io(self, output=None):
        """
        Convert the PRM object to an IO stream.

        Parameters
        ----------
        output : IO stream, optional
            The IO stream to write the PRM string to. Default is None.

        Returns
        -------
        str
            The PRM string.
        """
        return self.df.reset_index().astype(str).apply(lambda row: (' = ').join(row), axis=1)\
            .to_csv(output, header=None, index=None)

    def sel(self, *args):
        """
        Select specific PRM attributes and create a new PRM object.

        Parameters
        ----------
        *args : str
            The attribute names to select.

        Returns
        -------
        PRM
            The new PRM object with selected attributes.
        """
        return PRM(self.df.loc[[*args]])

    def __add__(self, other):
        """
        Add two PRM objects or a PRM object and a scalar.

        Parameters
        ----------
        other : PRM or scalar
            The PRM object or scalar to add.

        Returns
        -------
        PRM
            The resulting PRM object after addition.
        """
        import pandas as pd
        if isinstance(other, PRM):
            prm = pd.concat([self.df, other.df])
            # drop duplicates
            prm = prm.groupby(prm.index).last()
        else:
            prm = self.df + other
        return PRM(prm)

    def __sub__(self, other):
        """
        Subtract two PRM objects or a PRM object and a scalar.

        Parameters
        ----------
        other : PRM or scalar
            The PRM object or scalar to subtract.

        Returns
        -------
        PRM
            The resulting PRM object after subtraction.
        """
        import pandas as pd
        if isinstance(other, PRM):
            prm = pd.concat([self.df, other.df])
            # drop duplicates
            prm = prm.groupby(prm.index).last()
        else:
            prm = self.df - other
        return PRM(prm)

    def get(self, *args):
        """
        Get the values of specific PRM attributes.

        Parameters
        ----------
        *args : str
            The attribute names to get values for.

        Returns
        -------
        Union[Any, List[Any]]
            The values of the specified attributes. If only one attribute is requested, 
            return its value directly. If multiple attributes are requested, return a list of values.
        """
        out = [self.df.loc[[key]].iloc[0].values[0] for key in args]
        if len(out) == 1:
            return out[0]
        return out

    def shift_atime(self, lines, inplace=False):
        """
        Shift time in azimuth by a number of lines.

        Parameters
        ----------
        lines : float
            The number of lines to shift by.
        inplace : bool, optional
            Whether to modify the PRM object in-place. Default is False.

        Returns
        -------
        PRM
            The shifted PRM object or DataFrame. If 'inplace' is True, returns modified PRM object,
            otherwise, returns a new PRM with shifted times.
        """
        prm = self.sel('clock_start','clock_stop','SC_clock_start','SC_clock_stop') + lines/self.get('PRF')/86400.0
        if inplace:
            return self.set(prm)
        else:
            return prm

    # fitoffset.csh 3 3 freq_xcorr.dat
    # PRM.fitoffset(3, 3, offset_dat)
    # PRM.fitoffset(3, 3, matrix_fromfile='raw/offset.dat')
    @staticmethod
    def fitoffset(rank_rng, rank_azi, matrix=None, matrix_fromfile=None, SNR=20):
        """
        Estimates range and azimuth offsets for InSAR (Interferometric Synthetic Aperture Radar) data.

        Parameters
        ----------
        rank_rng : int
            Number of parameters to fit in the range direction.
        rank_azi : int
            Number of parameters to fit in the azimuth direction.
        matrix : numpy.ndarray, optional
            Array of range and azimuth offset estimates. Default is None.
        matrix_fromfile : str, optional
            Path to a file containing range and azimuth offset estimates. Default is None.
        SNR : int, optional
            Signal-to-noise ratio cutoff. Data points with SNR below this threshold are discarded.
            Default is 20.

        Returns
        -------
        prm : PRM object
            An instance of the PRM class with the calculated parameters.

        Raises
        ------
        Exception
            If both 'matrix' and 'matrix_fromfile' arguments are provided or if neither is provided.
        Exception
            If there are not enough data points to estimate the parameters.

        Usage
        -----
        The function estimates range and azimuth offsets for InSAR data based on the provided input.
        It performs robust fitting to obtain range and azimuth coefficients, calculates scale coefficients,
        and determines the range and azimuth shifts. The resulting parameters are then stored in a PRM object.

        Example
        -------
        fitoffset(3, 3, matrix_fromfile='raw/offset.dat')
        """
        import numpy as np
        import math

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

        rng_coef = PRM.robust_trend2d(rng, rank_rng)
        azi_coef = PRM.robust_trend2d(azi, rank_azi)

        # print MSE (optional)
        #rng_mse = PRM.robust_trend2d_mse(rng, rng_coef, rank_rng)
        #azi_mse = PRM.robust_trend2d_mse(azi, azi_coef, rank_azi)
        #print ('rng_mse_norm', rng_mse/len(rng), 'azi_mse_norm', azi_mse/len(azi))

        # range and azimuth data ranges
        scale_coef = [np.min(rng[:,0]), np.max(rng[:,0]), np.min(rng[:,1]), np.max(rng[:,1])]

        #print ('rng_coef', rng_coef)
        #print ('azi_coef', azi_coef)

        # now convert to range coefficients
        rshift = rng_coef[0] - rng_coef[1]*(scale_coef[1]+scale_coef[0])/(scale_coef[1]-scale_coef[0]) \
            - rng_coef[2]*(scale_coef[3]+scale_coef[2])/(scale_coef[3]-scale_coef[2])
        # now convert to azimuth coefficients
        ashift = azi_coef[0] - azi_coef[1]*(scale_coef[1]+scale_coef[0])/(scale_coef[1]-scale_coef[0]) \
            - azi_coef[2]*(scale_coef[3]+scale_coef[2])/(scale_coef[3]-scale_coef[2])
        #print ('rshift', rshift, 'ashift', ashift)

        # note: Python x % y expression and nympy results are different to C, use math function
        # use 'g' format for float values as in original GMTSAR codes to easy compare results
        prm = PRM().set(gformat=True,
                        rshift     =int(rshift) if rshift>=0 else int(rshift)-1,
                        sub_int_r  =math.fmod(rshift, 1)  if rshift>=0 else math.fmod(rshift, 1) + 1,
                        stretch_r  =rng_coef[1]*2/(scale_coef[1]-scale_coef[0]),
                        a_stretch_r=rng_coef[2]*2/(scale_coef[3]-scale_coef[2]),
                        ashift     =int(ashift) if ashift>=0 else int(ashift)-1,
                        sub_int_a  =math.fmod(ashift, 1)  if ashift>=0 else math.fmod(ashift, 1) + 1,
                        stretch_a  =azi_coef[1]*2/(scale_coef[1]-scale_coef[0]),
                        a_stretch_a=azi_coef[2]*2/(scale_coef[3]-scale_coef[2]),
                       )

        return prm

    def diff(self, other, gformat=True):
        """
        Compare the PRM object with another PRM object and return the differences.

        Parameters
        ----------
        other : PRM
            The other PRM object to compare with.
        gformat : bool, optional
            Whether to use 'g' format for float values. Default is True.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the differences between the two PRM objects.
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

    # note: only one dimension chunked due to sequential file reading 
    def read_SLC_int(self, intensity=False, chunksize=None):
        """
        Read SLC (Single Look Complex) data and compute the power of the signal.
        The method reads binary SLC data file, which contains alternating sequences of real and imaginary parts.
        It calculates the intensity of the signal and return it as a 2D numpy array.

        Returns
        -------
        xarray.DataArray
            2D array representing the power of the signal. The array shape corresponds to the dimensions of the SLC data.

        Notes
        -----
        This function uses a data factor (DFACT = 2.5e-07) from the GMTSAR code.
        The GMTSAR note indicates that the square of the intensity is used to match gips ihconv.
        The returned intensity data is flipped up-down ("shift data up if necessary") following the GMTSAR convention.

        Raises
        ------
        FileNotFoundError
            If the SLC file cannot be found.

        Example
        -------
        >>> import numpy as np
        >>> prm = PRM.from_file(filename)
        >>> amp = prm.read_SLC_int()
        """
        import xarray as xr
        import numpy as np
        import dask, dask.array
        import os

        # define lost class variables due to joblib
        if chunksize is None:
            chunksize = self.chunksize

        @dask.delayed
        def read_SLC_block(slc_filename, start, stop, intensity):
            # from GMTSAR code
            DFACT = 2.5e-07
            # Read a chunk of the SLC file
            # [real_0, imag_0, real_1, imag_1, real_2, imag_2, ...]
            #print ('    offset, shape', start*2, 2*(stop-start))
            # offset is measured in bytes
            data = np.memmap(slc_filename, dtype=np.int16, mode='r', offset=start*4, shape=(2*(stop-start),))
            #print ('        size', data.size)
            real_part = data[::2]
            imag_part = data[1::2]
            if not intensity:
                # return original complex data
                return (DFACT*(real_part + 1j * imag_part)).astype(np.complex64)
            # Calculate intensity (GMTSAR compatible while it names intensity as amplitude)
            #return (DFACT*np.abs(real_part + 1j * imag_part))**2
            return ((DFACT*real_part)**2 + (DFACT*imag_part)**2).astype(np.float32)

        prm = PRM.from_file(self.filename)
        # num_patches multiplier is omitted
        slc_filename, xdim, ydim = prm.get('SLC_file', 'num_rng_bins', 'num_valid_az')
        dirname = os.path.dirname(self.filename)
        slc_filename = os.path.join(dirname, slc_filename)
        #print (slc_filename, ydim, xdim)

        blocksize = chunksize*xdim
        blocks = int(np.ceil(ydim * xdim / blocksize))
        #print ('chunks', chunks, 'chunksize', chunksize)
        # Create a lazy Dask array that reads chunks of the SLC file
        lazy_arrays = []
        for i in range(blocks):
            start = i * blocksize
            stop = min((i+1) * blocksize, ydim * xdim)
            #print ('start, stop, shape', start, stop, (stop-start))
            # use proper output data type for complex data and intensity
            block = dask.array.from_delayed(read_SLC_block(slc_filename, start, stop, intensity), shape=((stop-start),),
                dtype=np.float32 if intensity else np.complex64)
            lazy_arrays.append(block)
            del block
        # Concatenate the chunks together
        # Data file can include additional data outside of the specified dimensions
        data = dask.array.concatenate(lazy_arrays).reshape((-1, xdim))[:ydim,:]
        del lazy_arrays
        #amp = xr.DataArray(dask.array.flipud(amp), coords={'y': np.arange(ydim) + 0.5, 'x': np.arange(xdim) + 0.5})
        data = xr.DataArray(data, coords={'y': np.arange(ydim) + 0.5, 'x': np.arange(xdim) + 0.5})
        # set the short name instead of autogenerated one
        return data.rename('SLC')

    @staticmethod
    def goldstein_filter_parallel(data, corr, psize):
        import xarray as xr
        import numpy as np
        import dask
    
        def apply_pspec(data, alpha):
            # NaN is allowed value
            assert not(alpha < 0), f'Invalid parameter value {alpha} < 0'
            wgt = np.power(np.abs(data)**2, alpha / 2)
            data = wgt * data
            return data

        def make_wgt(nxp, nyp):
            # Create arrays of horizontal and vertical weights
            wx = 1.0 - np.abs(np.arange(nxp // 2) - (nxp / 2.0 - 1.0)) / (nxp / 2.0 - 1.0)
            wy = 1.0 - np.abs(np.arange(nyp // 2) - (nyp / 2.0 - 1.0)) / (nyp / 2.0 - 1.0)
            # Compute the outer product of wx and wy to create the top-left quadrant of the weight matrix
            quadrant = np.outer(wy, wx)
            # Create a full weight matrix by mirroring the quadrant along both axes
            wgt = np.block([[quadrant, np.flip(quadrant, axis=1)],
                            [np.flip(quadrant, axis=0), np.flip(np.flip(quadrant, axis=0), axis=1)]])
            return wgt

        def goldstein_filter(data, corr, wgt, psize):
            """
            Apply the Goldstein adaptive filter to the given data.

            Args:
                data: 2D numpy array of complex values representing the data to be filtered.
                corr: 2D numpy array of correlation values. Must have the same shape as `data`.

            Returns:
                2D numpy array of filtered data.
            """
            # Calculate alpha
            alpha = 1 - (wgt * corr).sum() / wgt.sum()
            data = np.fft.fft2(data, s=(psize,psize))
            data = apply_pspec(data, alpha)
            data = np.fft.ifft2(data, s=(psize,psize))
            return wgt * data

        def apply_goldstein_filter(data, corr, psize):
            # Create an empty array for the output
            out = np.zeros(data.shape, dtype=np.complex64)
            # ignore processing for empty chunks 
            if np.all(np.isnan(data)):
                return out
            # Create the weight matrix
            wgt_matrix = make_wgt(psize, psize)
            # Iterate over windows of the data
            for i in range(0, data.shape[0] - psize, psize // 2):
                for j in range(0, data.shape[1] - psize, psize // 2):
                    # Create proocessing windows
                    data_window = data[i:i+psize, j:j+psize]
                    corr_window = corr[i:i+psize, j:j+psize]
                    wgt_window = wgt_matrix[:data_window.shape[0],:data_window.shape[1]]
                    # Apply the filter to the window
                    filtered_window = goldstein_filter(data_window, corr_window, wgt_window, psize)
                    # Add the result to the output array
                    slice_i = slice(i, min(i + psize, out.shape[0]))
                    slice_j = slice(j, min(j + psize, out.shape[1]))
                    out[slice_i, slice_j] += filtered_window[:slice_i.stop - slice_i.start, :slice_j.stop - slice_j.start]
            return out

        # Apply function with overlap; psize//2 overlap is not enough (some empty lines produced)
        # use complex data and real correlation
        phase_filtered = dask.array.map_overlap(lambda data, corr: apply_goldstein_filter(data, corr, psize),
                                       data.data,
                                       corr.data,
                                       depth=psize//2 + 2,
                                       dtype=np.float32, 
                                       meta=np.array(()))
        # Calculate the phase and fix chunksizes
        return xr.DataArray(np.arctan2(phase_filtered.imag, phase_filtered.real),coords=corr.coords)\
                .chunk(corr.chunksizes).rename('phase')

    @staticmethod
    def correlation(A1, A2, amp):
        import xarray as xr
        import numpy as np
        # constant from GMTSAR code
        thresh = 5.e-21
        a = A1 * A2
        corr = xr.where(a > 0, amp / np.sqrt(a), 0)
        corr = xr.where(corr < 0, 0, corr)
        corr = xr.where(corr > 1, 1, corr)
        # mask too low amplitude areas as invalid
        # amp1 and amp2 chunks are high for SLC, amp has normal chunks for NetCDF
        return xr.where(a >= thresh, corr, np.nan).chunk(a.chunksizes).rename('phase')

    # see about correlation filter
    # https://github.com/gmtsar/gmtsar/issues/86
    def intf(self, other, basedir, topo_ra_fromfile, basename=None, wavelength=200, psize=32, \
            coarsen=(1,4), func=None, weight=None, chunksize=None, debug=False):
        """
        Perform interferometric processing on the input SAR data.

        Parameters
        ----------
        other : PRM
            Another instance of the PRM class representing the second SAR data.
        basedir : str
            Base directory for storing the processed data files.
        topo_ra_fromfile : str
            Path to the topo_ra file used for the calculation.
        basename : str, optional
            Base name for the output files. If not provided, a default name will be generated.
        wavelength : int, optional
            Wavelength of the SAR data in meters. Default is 200 m.
        psize : int, optional
            Werner/Goldstein filter window size in pixels. Default is 32 pixels.
        func : callable, optional
            Custom function to apply on the processed data arrays. Default is None.
        chunksize : int or dict, optional
            Chunk size for dask arrays. Default is None.
        debug : bool, optional
            Enable debug mode. Default is False.

        Returns
        -------
        None

        Notes
        -----
        This method performs interferometric processing on the input SAR data using GMTSAR tools.
        It generates various intermediate and final files, including amplitude and phase data,
        and applies filtering and convolution operations to produce the desired output.
        """
        import os
        import numpy as np
        import xarray as xr
        import dask.array

        if not isinstance(other, PRM):
            raise Exception('Argument "other" should be PRM class instance')

        if coarsen is None and psize is not None:
            raise Exception('Argument "coarsen" can be None only when "psize" is None too')

        # define lost class variables due to joblib
        if chunksize is None:
            chunksize = self.chunksize

        # define basename from two PRM names
        subswath = os.path.basename(self.filename)[16:18]
        if basename is None:
            # S1_20171111_ALL_F3.PRM -> F3
            date1    = os.path.basename(self.filename)[3:11]
            date2    = os.path.basename(other.filename)[3:11]
            basename = f'{subswath}_{date1}_{date2}_'
            #print ('basename', basename)

        # make full file name, use workaround for 'weight' argument name defined without extension
        fullname = lambda name: os.path.join(basedir, basename + name if name[-4:]=='.grd' else f'{subswath}_{name}.grd')

        # weight can be None, xarray dataarray, NetCDF file name
        if weight is not None and isinstance(weight, str):
            if debug:
                print (f'DEBUG: intf weight from file {weight}')
            weight = xr.open_dataarray(fullname(weight), engine=self.engine, chunks=chunksize)
            # revert fake dimension names 
            if weight.dims == ('a', 'r'):
                weight = weight.rename({'a': 'y', 'r': 'x'})

        # prepare PRMs for the calculation below
        other.set(self.SAT_baseline(other, tail=9))
        self.set(self.SAT_baseline(self).sel('SC_height','SC_height_start','SC_height_end'))

        # for topo_ra use relative path from PRM files directory
        # use imag.grd=bf for GMT native, C-binary format
        # os.path.join(basedir, f'{subswath}_topo_ra.grd')
        if os.path.exists(fullname('real.grd')):
            os.remove(fullname('real.grd'))
        if os.path.exists(fullname('imag.grd')):
            os.remove(fullname('imag.grd'))
        self.phasediff(other, topo_ra_fromfile=topo_ra_fromfile,
                       imag_tofile=fullname('imag.grd'),
                       real_tofile=fullname('real.grd'),
                       debug=debug)

        # original SLC (do not flip vertically)
        amp1 = self.read_SLC_int(intensity=True)
        amp2 = other.read_SLC_int(intensity=True)
        # phasediff tool output files (flip vertically)
        imag = xr.open_dataarray(fullname('imag.grd'), engine=self.engine, chunks=chunksize)
        imag.data = dask.array.flipud(imag)
        real = xr.open_dataarray(fullname('real.grd'), engine=self.engine, chunks=chunksize)
        real.data = dask.array.flipud(real)
        #real.shape, imag.shape (5484, 21572) (5484, 21572)
        #print ('DEBUG X1 real.shape, imag.shape', real.shape, imag.shape)

        # anti-aliasing filter for multi-looking, wavelength can be None
        imag = self.antialiasing_downscale(imag, weight=weight, wavelength=wavelength, coarsen=coarsen, debug=debug)
        real = self.antialiasing_downscale(real, weight=weight, wavelength=wavelength, coarsen=coarsen, debug=debug)
        amp1 = self.antialiasing_downscale(amp1, weight=weight, wavelength=wavelength, coarsen=coarsen, debug=debug)
        amp2 = self.antialiasing_downscale(amp2, weight=weight, wavelength=wavelength, coarsen=coarsen, debug=debug)

        # calculate amplitude of interferogram
        amp = np.sqrt(real**2 + imag**2)
        # calculate masked correlation
        corr = self.correlation(amp1, amp2, amp)
        # cleanup
        del amp1, amp2, amp
        #chunksize=(512, 5393)
        #print ('DEBUG X2 corr', corr)

        # Apply Goldstein filter function after multi-looking
        if debug:
            print ('DEBUG: intf apply Goldstein filter with size', psize is not None, psize)
        if psize is not None:
            phase = self.goldstein_filter_parallel((real + 1j * imag), corr, psize=psize)
        else:
            phase = np.arctan2(imag, real)
        # cleanup
        del real, imag
        #chunksize=(512, 5393)
        #print ('DEBUG X3 phase', phase)

        if func is not None:
            corr = func(corr)
        if os.path.exists(fullname('corr.grd')):
            os.remove(fullname('corr.grd'))
        encoding = {'corr': self.compression(corr.shape, chunksize=chunksize)}
        #print ('DEBUG X2', corr_da)
        # rename to save lazy NetCDF preventing broken coordinates (y,y) 
        corr.rename('corr').rename({'y': 'a', 'x': 'r'}).to_netcdf(fullname('corr.grd'), encoding=encoding, engine=self.engine)

        if func is not None:
            phase = func(phase)
        if os.path.exists(fullname('phasefilt.grd')):
            os.remove(fullname('phasefilt.grd'))
        encoding = {'phase': self.compression(phase.shape, chunksize=chunksize)}
        #print ('DEBUG X3', phase_da)
        # mask phase using masked correlation
        # rename to save lazy NetCDF preventing broken coordinates (y,y)
        phase.rename('phase').rename({'y': 'a', 'x': 'r'}).to_netcdf(fullname('phasefilt.grd'), encoding=encoding, engine=self.engine)

        # cleanup
        del corr, phase
        for name in ['real.grd', 'imag.grd']:
            filename = fullname(name)
            if os.path.exists(filename):
                os.remove(filename)

        #print ('DEBUG X4')
        return

    # see make_gaussian_filter.c for the original code
    def pixel_size(self, grid=1):
        """
        Calculate azimuth and range pixel size in meters.

        Parameters
        ----------
        grid : tuple or xarray.DataArray, optional
            A pair of x, y grid decimation coefficients or a 2D or 3D Xarray DataArray representing the decimation grid.
            The default is 1 for the input Sentinel-1 SLC grid.

        Returns
        -------
        tuple
            A tuple containing the azimuth and range pixel sizes in meters.

        Notes
        -----
        This method calculates the pixel sizes in the azimuth and range directions based on the provided parameters.
        The azimuth pixel size is determined using the spacecraft velocity, height, and pulse repetition frequency (PRF),
        while the range pixel size is calculated based on the speed of light and the range sampling rate.
        """
        import xarray as xr
        import numpy as np
        from scipy.constants import speed_of_light

        # define grid spacing
        if isinstance(grid, (xr.DataArray, xr.Dataset)):
            assert self.is_ra(grid), 'Raster needs to be defined in radar coordinates'
            dy = grid.y.diff('y')[0].item()
            dx = grid.x.diff('x')[0].item()
        elif isinstance(grid, (list, tuple)):
            dy, dx = grid
        else:
            dy = dx = grid
            
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
        return (dy * azi_px_size, dx * rng_px_size)

    # TODO: use PRM parameters to define config parameters
    def snaphu_config(self, defomax=0, **kwargs):
        """
        Generate a configuration file for Snaphu phase unwrapping.

        Parameters
        ----------
        defomax : int, optional
            Maximum number of deformation cycles allowed. Default is 0.
        **kwargs
            Additional configuration parameters in key-value format.

        Returns
        -------
        str
            The generated configuration file as a string.

        Examples
        --------
        Get default SNAPHU config:
        config = sbas.snaphu_config()

        Get custom SNAPHU config with added tiling options:
        config = sbas.snaphu_config(defomax=DEFOMAX, NTILEROW=1, NTILECOL=2, ROWOVRLP=200, COLOVRLP=200)

        Notes
        -----
        This method generates a configuration file for Snaphu, a phase unwrapping software.
        The configuration file includes basic parameters and allows customization by passing additional parameters.

        Custom Parameters
        -----------------
        Additional configuration parameters can be passed as keyword arguments (e.g., `parameter_name=value`).
        Boolean values should be provided as `True` or `False`.
        """
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

