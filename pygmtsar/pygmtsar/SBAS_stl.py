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
from .tqdm_dask import tqdm_dask

class SBAS_stl(SBAS_incidence):

    @staticmethod
    def stl(ts, dt, dt_periodic, periods, incremental=True, robust=False):
        """
        Perform Seasonal-Trend decomposition using LOESS (STL) on the input time series data.

        The function performs the following steps:
        1. Check for NaN values in the input time series and return arrays filled with NaNs if found.
        2. If the input time series data is incremental, calculate the real values first.
        3. Create an interpolation function using the input time series and corresponding time values.
        4. Interpolate the time series data for the periodic time values.
        5. Perform STL decomposition using the interpolated time series data.
        6. Return the trend, seasonal, and residual components of the decomposed time series.

        Parameters
        ----------
        ts : numpy.ndarray
            Input time series data.
        dt : numpy.ndarray
            Corresponding time values for the input time series data.
        dt_periodic : numpy.ndarray
            Periodic time values for interpolation.
        periods : int
            Number of periods for seasonal decomposition.
        incremental : bool, optional
            Whether the input time series data is incremental (default is True).
        robust : bool, optional
            Whether to use a robust fitting procedure for the STL decomposition (default is False).

        Returns
        -------
        numpy.ndarray
            Trend component of the decomposed time series.
        numpy.ndarray
            Seasonal component of the decomposed time series.
        numpy.ndarray
            Residual component of the decomposed time series.

        """
        import numpy as np
        from scipy.interpolate import interp1d
        from statsmodels.tsa.seasonal import STL

        # Check for NaNs in the input time series; this check is faster than letting STL handle NaNs
        if np.any(np.isnan(ts)):
            nodata = np.nan * np.zeros(dt_periodic.size)
            return nodata, nodata, nodata

        # Calculate the real values from the incremental input time series
        if incremental:
            ts = np.cumsum(ts)

        # Create an interpolation function for the input time series and time values
        interp_func = interp1d(dt, ts, kind='nearest', fill_value='extrapolate', assume_sorted=True)
        # Interpolate the time series data for the periodic time values
        ts = interp_func(dt_periodic)

        # Perform STL decomposition on the interpolated time series data
        stl = STL(ts, period=periods, robust=robust)
        res = stl.fit()

        return res.trend, res.seasonal, res.resid

    # Aggregate data for varying frequencies (e.g., 12+ days for 6 days S1AB images interval)
    # Use frequency strings like '1W' for 1 week, '2W' for 2 weeks, '10d' for 10 days, '1M' for 1 month, etc.
    # Calculate cumulative sum along the date dimension when incremental=True
    def stl_parallel(self, dates, data='disp', freq='W', periods=52, incremental=True, robust=False, chunksize=None, interactive=False):
        """
        Perform Seasonal-Trend decomposition using LOESS (STL) on the input time series data in parallel.

        The function performs the following steps:
        1. Convert the 'date' coordinate to valid dates.
        2. Unify date intervals to a specified frequency (e.g., weekly) for a mix of time intervals.
        3. Apply the sbas_stl function in parallel using xarray's apply_ufunc and Dask.
        4. Rename the output date dimension to match the original irregular date dimension.
        5. Return the STL decomposition results as an xarray Dataset.

        Parameters
        ----------
        self : SBAS
            Instance of the SBAS class.
        dates : numpy.ndarray
            Array of datetime64 values corresponding to the input time series data.
        data : str or xarray.DataArray, optional
            Input time series data as an xarray DataArray or the name of the grid (default is 'disp').
        freq : str, optional
            Frequency string for unifying date intervals (default is 'W' for weekly).
        periods : int, optional
            Number of periods for seasonal decomposition (default is 52).
        incremental : bool, optional
            Whether to calculate data cumsum for the date dimension when True (default is True).
        robust : bool, optional
            Whether to use a slower robust fitting procedure for the STL decomposition (default is False).
        chunksize : int, optional
            The chunk size to be used for dimensions (default is None, which will use the class-level chunksize).
        interactive : bool, optional
            Whether to return the model interactively or save it to NetCDF file (default is False).

        Returns
        -------
        xarray.Dataset or None
            An xarray Dataset containing the trend, seasonal, and residual components of the decomposed time series,
            or None if the results are saved to a file.

        See Also
        --------
        statsmodels.tsa.seasonal.STL : Seasonal-Trend decomposition using LOESS
            https://www.statsmodels.org/dev/generated/statsmodels.tsa.seasonal.STL.html
        """
        import xarray as xr
        import numpy as np
        import pandas as pd
        import dask
        import os

        if chunksize is None:
            # smaller chunks are better for large 3D grids processing
            chunksize = self.chunksize

        if isinstance(data, str):
            data = self.open_grids(dates, data, func=[self.vertical_displacement_mm], geocode=True, chunksize=chunksize)
        elif isinstance(data, xr.DataArray):
            pass
        else:
            raise Exception('Invalid input: The "data" parameter should be of type xarray.DataArray or a string representing a grid name.')
        
        # convert coordinate to valid dates
        data['date'] = pd.to_datetime(data.date)
    
        # Unify date intervals; using weekly intervals should be suitable for a mix of 6 and 12 days intervals
        dates_weekly = pd.date_range(dates[0], dates[-1], freq=freq)
        dt_weekly = xr.DataArray(dates_weekly, dims=['date'])
        # The following line of code is not efficient because the "date" dimension is changed
        # Note: The bottleneck library should be installed for this line of code to work
        # cumdisp_weekly = data.cumsum('date').interp(date=dates_weekly, method='nearest', assume_sorted=True)

        # Output "date" is different from the input "date" dimension
        trend, seasonal, resid = xr.apply_ufunc(
            self.stl,
            data.chunk(dict(date=-1)),
            data.date.astype(np.int64),
            input_core_dims=[['date'],['date']],
            output_core_dims=[['date2'], ['date2'], ['date2']],
            dask='parallelized',
            vectorize=True,
            output_dtypes=[np.float32, np.float32, np.float32],
            output_sizes={'date2': len(dates_weekly)},
            kwargs={'dt_periodic': dt_weekly.astype(np.int64),
                    'periods': periods,
                    'incremental': incremental,
                    'robust': robust}
        )
        # rename regular date dimension to "date" as original irregular date
        model = xr.Dataset({'trend': trend, 'seasonal': seasonal, 'resid': resid},
                           coords={'date2': dt_weekly.values}).rename({'date2': 'date'})
        if interactive:
            return model

        # Define the encoding and choose the appropriate storage scheme
        model_filename = self.get_filenames(None, None, f'stl')
        if os.path.exists(model_filename):
            os.remove(model_filename)
        netcdf_compression = self.compression(chunksize=(1, chunksize, chunksize))
        delayed = model.to_netcdf(model_filename,
                     engine=self.engine,
                     encoding={varname: netcdf_compression for varname in model.data_vars},
                     compute=False)
        tqdm_dask(dask.persist(delayed), desc='Seasonal-Trend decomposition using LOESS')
