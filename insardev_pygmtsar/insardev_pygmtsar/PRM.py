# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from insardev_toolkit import datagrid
from .PRM_gmtsar import PRM_gmtsar

class PRM(datagrid, PRM_gmtsar):
    
    int_types = ['num_valid_az', 'num_rng_bins', 'num_patches', 'bytes_per_line', 'good_bytes_per_line', 'num_lines','SC_identity']

    @staticmethod
    def SC_timestamp(SC_clock):
        from datetime import datetime, timedelta

        # extract year and Julian day with fractional part
        year = int(SC_clock // 1000)
        julian_day = SC_clock % 1000  # Keep fractional part of the day

        # split integer Julian day and fractional day
        integer_julian_day = int(julian_day)
        fractional_day = julian_day - integer_julian_day

        # convert integer and fraction parts of Julian day to datetime
        timestamp = datetime(year, 1, 1) + timedelta(days=integer_julian_day) + timedelta(days=fractional_day)

        return timestamp

    @staticmethod
    def to_numeric_or_original(val):
        if isinstance(val, str):
            try:
                float_val = float(val)
                return float_val
            except ValueError:
                return val
        return val

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
        if isinstance(prm_string, bytes):
            # for cases like
            #return PRM.from_str(os.read(pipe2[0],int(10e6))
            prm_string = prm_string.decode('utf-8')
        # for cases like
        # return PRM.from_str(os.read(pipe2[0],int(10e6).decode('utf8'))
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

        df = pd.read_csv(prm, sep=r'\s+=\s+', header=None, names=['name', 'value'], engine='python').set_index('name')
        df['value'] = df['value'].map(PRM.to_numeric_or_original)

        return PRM(df)

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

        # Initialize an empty DataFrame if prm is None
        if prm is None:
            _prm = pd.DataFrame(columns=['name', 'value'])
        elif isinstance(prm, pd.DataFrame):
            _prm = prm.reset_index()
        else:
            _prm = prm.df.reset_index()

        # Convert values to numeric where possible, keep original value otherwise
        _prm['value'] = _prm['value'].map(PRM.to_numeric_or_original)

        # Set the DataFrame for the PRM object
        self.df = _prm[['name', 'value']].drop_duplicates(keep='last').set_index('name')
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
    def set(self, prm=None, **kwargs):
        """
        Set PRM values.

        Parameters
        ----------
        prm : PRM, optional
            The PRM object to set values from. Default is None.
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
            self.df.loc[key] = value
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
        # update internal filename after saving with the new filename
        self.filename = prm
        return self

    def update(self, debug=False):
        """
        Save PRM file to disk.

        Parameters
        ----------
        debug : bool, optional
            Whether to enable debug mode. Default is False.

        Returns
        -------
        PRM
            The updated PRM object.
        """
        if self.filename is None:
            raise Exception('PRM is not created from file, use to_file() method instead')

        if debug:
            print ('DEBUG:', self)

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
        vals = [self.df.loc[[key]].iloc[0].values[0] for key in args]
        out = [int(v) if k in self.int_types else v for (k, v) in zip(args,vals)]
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

    def fix_aligned(self):
        """
        Correction for the range and azimuth shifts of the re-aligned SLC images (fix_prm_params() in GMTSAR)
        """
        from scipy import constants
        # constant from GMTSAR code
        #SOL = 299792456.0
    
        delr = constants.speed_of_light / self.get('rng_samp_rate') / 2
        #delr = SOL / self.get('rng_samp_rate') / 2
        near_range = self.get('near_range') + \
            (self.get('st_rng_bin') - self.get('chirp_ext') + self.get('rshift') + self.get('sub_int_r') - 1)* delr
    
        SC_clock_start = self.get('SC_clock_start') + \
            (self.get('ashift') + self.get('sub_int_a')) / (self.get('PRF') * 86400.0) + \
            (self.get('nrows') - self.get('num_valid_az')) / (2 * self.get('PRF') * 86400.0)
    
        SC_clock_stop = SC_clock_start + \
            (self.get('num_valid_az') * self.get('num_patches')) / (self.get('PRF') * 86400.0)
    
        return self.set(near_range=near_range,
                        SC_clock_start=SC_clock_start,
                        SC_clock_stop=SC_clock_stop)

    # note: only one dimension chunked due to sequential file reading 
    def read_SLC_int(self, scale=2.5e-07, shape=None):
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
        import dask
        import dask.array as da
        import os
        import warnings

        @dask.delayed
        def read_SLC_block(slc_filename, start, stop):
            # Read a chunk of the SLC file
            # [real_0, imag_0, real_1, imag_1, real_2, imag_2, ...]
            #print ('    offset, shape', start*2, 2*(stop-start))
            # offset is measured in bytes
            return np.memmap(slc_filename, dtype=np.int16, mode='r', offset=start*4, shape=(2*(stop-start),))

        prm = PRM.from_file(self.filename)
        # num_patches multiplier is omitted
        #slc_filename, xdim, ydim = prm.get('SLC_file', 'num_rng_bins', 'num_valid_az')
        slc_filename, xdim, ydim, rshift, ashift = prm.get('SLC_file', 'num_rng_bins', 'num_valid_az', 'rshift', 'ashift')
                
        dirname = os.path.dirname(self.filename)
        slc_filename = os.path.join(dirname, slc_filename)
        #print (slc_filename, ydim, xdim)

        blocksize = self.netcdf_chunksize * xdim
        blocks = int(np.ceil(ydim * xdim / blocksize))
        #print ('chunks', chunks, 'chunksize', chunksize)
        # Create a lazy Dask array that reads chunks of the SLC file
        res = []
        ims = []
        for i in range(blocks):
            start = i * blocksize
            stop = min((i+1) * blocksize, ydim * xdim)
            #assert 0, f'SLC: start={start}, stop={stop}'
            #print ('start, stop, shape', start, stop, (stop-start))
            # use proper output data type for complex data and intensity
            block = da.from_delayed(read_SLC_block(slc_filename, start, stop),
                shape=(2*(stop-start),), dtype=np.int16)
            res.append(block[::2])
            ims.append(block[1::2])
            del block
        # Concatenate the chunks together
        # Data file can include additional data outside of the specified dimensions
        re = da.concatenate(res).reshape((-1, xdim))[:ydim,:]
        im = da.concatenate(ims).reshape((-1, xdim))[:ydim,:]
        del res, ims
        
        assert re.shape == (ydim, xdim), f'Originated re shape ({ydim},{xdim}), but got {re.shape}'
        assert im.shape == (ydim, xdim), f'Originated im width ({ydim},{xdim}), but got {im.shape}'

        # pad to the specified reference frame
        if shape is not None:
            if shape[1] is not None and shape[1] > xdim:
                rpad = da.zeros((ydim, shape[1] - xdim), dtype=np.int16)
                re = da.concatenate([re, rpad], axis=1)
                im = da.concatenate([im, rpad], axis=1)
            elif shape[1] is not None and shape[1] < xdim:
                re = re[:,:shape[1]]
                im = im[:,:shape[1]]
            if shape[0] is not None and shape[0] > ydim:
                apad = da.zeros((shape[0] - ydim, shape[1] if shape[1] is not None else xdim), dtype=np.int16)
                re = da.concatenate([re, apad], axis=0)
                im = da.concatenate([im, apad], axis=0)
            elif shape[0] is not None and shape[0] < ydim:
                re = re[:shape[0]]
                im = im[:shape[0]]

        coords = {'a': np.arange(ydim if shape is None or shape[0] is None else shape[0]) + 0.5,
                  'r': np.arange(xdim if shape is None or shape[1] is None else shape[1]) + 0.5}

        #coords = {'y': np.arange(ydim) + 0.5, 'x': np.arange(xdim) + 0.5}
        # perform aligning on-the-fly using coordinates
        #coords = {'y': np.arange(ydim) + ashift + 0.5, 'x': np.arange(xdim) + rshift + 0.5}
        re = xr.DataArray(re, coords=coords).rename('re')
        im = xr.DataArray(im, coords=coords).rename('im')
        #if intensity:
        #    return ((scale or 1)*re.astype(np.float32))**2 + ((scale or 1)*im.astype(np.float32))**2
        if scale is not None:
            return scale * (xr.merge([re, im]).astype(np.float32))
        return xr.merge([re, im])

    def read_LED(self):
        """
        Read an associated LED file and extract the metadata and data into a dictionary and a DataFrame.
    
        Parameters:
        None
    
        Returns:
        tuple: A tuple containing the metadata dictionary and the data DataFrame.
    
        Metadata Dictionary:
        - 'nd' (str): Number of given points in the list.
        - 'iy' (str): Year.
        - 'id' (str): Julian day.
        - 'isec' (str): Seconds of the day.
        - 'idsec' (str): Delta seconds.
    
        Data DataFrame:
        - 'iy' (int): Year.
        - 'id' (int): Julian day.
        - 'isec' (float): Seconds of the day.
        - 'px' (float): X-coordinate.
        - 'py' (float): Y-coordinate.
        - 'pz' (float): Z-coordinate.
        - 'vx' (float): X-velocity.
        - 'vy' (float): Y-velocity.
        - 'vz' (float): Z-velocity.
        - 'time' (datetime): Combined datetime.
        """
        import pandas as pd
        from datetime import datetime, timedelta
    
        led_file = self.get('led_file')
        with open(led_file, 'r') as file:
            first_line = file.readline()
        columns = ['nd', 'iy', 'id', 'isec', 'idsec']
        meta = dict(zip(columns, first_line.split()))
    
        # Define the column headers based on the structure in the code
        headers = ['iy', 'id', 'isec', 'px', 'py', 'pz', 'vx', 'vy', 'vz']
        # Read the file with pandas, specifying space as the delimiter and skipping the first row
        df = pd.read_csv(led_file, delimiter=' ', header=None, skiprows=1, names=headers, index_col=False)
        df.attrs = meta
        
        # add a time column
        #df['time'] = df.apply(lambda row: datetime(int(row['iy']), 1, 1) + timedelta(days=int(row['id']), seconds=float(row['isec'])), axis=1)
        # add seconds from Jan 1
        df['clock'] = (24 * 60 * 60) * df['id'] + df['isec']

        return df

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
        prm = PRM().set(rshift     =int(rshift) if rshift>=0 else int(rshift)-1,
                        sub_int_r  =math.fmod(rshift, 1)  if rshift>=0 else math.fmod(rshift, 1) + 1,
                        stretch_r  =rng_coef[1]*2/(scale_coef[1]-scale_coef[0]),
                        a_stretch_r=rng_coef[2]*2/(scale_coef[3]-scale_coef[2]),
                        ashift     =int(ashift) if ashift>=0 else int(ashift)-1,
                        sub_int_a  =math.fmod(ashift, 1)  if ashift>=0 else math.fmod(ashift, 1) + 1,
                        stretch_a  =azi_coef[1]*2/(scale_coef[1]-scale_coef[0]),
                        a_stretch_a=azi_coef[2]*2/(scale_coef[3]-scale_coef[2]),
                       )

        return prm
