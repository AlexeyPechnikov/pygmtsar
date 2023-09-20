# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_tidal import Stack_tidal
from .tqdm_dask import tqdm_dask

class Stack_stl(Stack_tidal):

    def stack_vs_cube(self, data_pairs, data):
        """
        Calculate difference between pairs and dates
    
        Use to calculate solution vs pair unwrapped phases difference as
        diff = sbas.stack_vs_cube(phase_unwrap, solution)
        error = np.sqrt(sbas.wrap(diff)**2).sum('pair')
        """
        import xarray as xr

        # extract pairs
        pairs = self.get_pairs(data_pairs)
        # calculate differences between end and start dates for all the pairs
        error_pairs = []
        for rec in pairs.itertuples():
            error_pair = data.sel(date=rec.rep) - data.sel(date=rec.ref) - data_pairs.isel(pair=rec.Index)
            error_pairs.append(error_pair)
        # form 3D stack
        return xr.concat(error_pairs, dim='pair').assign_coords({'pair': data_pairs.pair})

    def cube_trend(self, data, deg=1):
        import xarray as xr

        trend = xr.polyval(data.date, data.polyfit('date', deg).polyfit_coefficients)
        return trend
    
    @staticmethod
    def stl1d(ts, dt, dt_periodic, periods, robust=False):
        """
        Perform Seasonal-Trend decomposition using LOESS (STL) on the input time series data.

        The function performs the following steps:
        1. Check for NaN values in the input time series and return arrays filled with NaNs if found.
        2. Create an interpolation function using the input time series and corresponding time values.
        3. Interpolate the time series data for the periodic time values.
        4. Perform STL decomposition using the interpolated time series data.
        5. Return the trend, seasonal, and residual components of the decomposed time series.

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
    def cube_stl(self, data, freq='W', periods=52, robust=False, chunksize=None, interactive=False):
        """
        Perform Seasonal-Trend decomposition using LOESS (STL) on the input time series data in parallel.

        The function performs the following steps:
        1. Convert the 'date' coordinate to valid dates.
        2. Unify date intervals to a specified frequency (e.g., weekly) for a mix of time intervals.
        3. Apply the Stack.stl1d function in parallel using Dask.
        4. Rename the output date dimension to match the original irregular date dimension.
        5. Return the STL decomposition results as an xarray Dataset.

        Parameters
        ----------
        self : Stack
            Instance of the Stack class.
        dates : numpy.ndarray
            Array of datetime64 values corresponding to the input time series data.
        data : xarray.DataArray
            Input time series data as an xarray DataArray.
        freq : str, optional
            Frequency string for unifying date intervals (default is 'W' for weekly).
        periods : int, optional
            Number of periods for seasonal decomposition (default is 52).
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

        Examples
        --------
        Use on (date,lat,lon) and (date,y,x) grids to return the results or store them on disk:
        stack.cube_stl(disp, interactive=True)
        stack.cube_stl(disp)

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
            # see lstsq() for details
            chunksize = self.chunksize // 4

        if isinstance(data, xr.DataArray):
            pass
        else:
            raise Exception('Invalid input: The "data" parameter should be of type xarray.DataArray.')
    
        # convert coordinate to valid dates
        dates = pd.to_datetime(data.date)

        dim0, dim1, dim2 = data.dims
        assert dim0 == 'date', 'The first data dimension should be date'
    
        # original dates
        dt = pd.to_datetime(data.date).astype(np.int64)
        # Unify date intervals; using weekly intervals should be suitable for a mix of 6 and 12 days intervals
        dates_weekly = pd.date_range(dates[0], dates[-1], freq=freq)
        dt_weekly = xr.DataArray(dates_weekly, dims=['date'])
        dt_periodic = dt_weekly.astype(np.int64)
        # The following line of code is not efficient because the "date" dimension is changed
        # Note: The bottleneck library should be installed for this line of code to work
        # disp_weekly = data.interp(date=dates_weekly, method='nearest', assume_sorted=True)

        def stl_block(lats, lons):
            # use external variables dt, dt_periodic, periods, robust
            # 3D array
            data_block = data.isel({dim1: lats, dim2: lons}).chunk(-1).compute(n_workers=1).values.transpose(1,2,0)
            # Vectorize vec_lstsq
            #vec_stl = np.vectorize(lambda data: self.stl(data, dt, dt_periodic, periods, robust), signature='(n)->(3,m)')
            vec_stl = np.vectorize(lambda data: self.stl1d(data, dt, dt_periodic, periods, robust), signature='(n)->(m),(m),(m)')
            # Apply vec_lstsq to data_block and revert the original dimensions order
            block = vec_stl(data_block)
            del vec_stl, data_block
            return np.asarray(block).transpose(0,3,1,2)

        # split to square chunks
        # use indices instead of the coordinate values to prevent the weird error raising occasionally:
        # "Reindexing only valid with uniquely valued Index objects"
        lats_blocks = np.array_split(np.arange(data[dim1].size), np.arange(0, data[dim1].size, chunksize)[1:])
        lons_blocks = np.array_split(np.arange(data[dim2].size), np.arange(0, data[dim2].size, chunksize)[1:])

        blocks_total = []
        for lats_block in lats_blocks:
            blocks = []
            for lons_block in lons_blocks:
                block = dask.array.from_delayed(dask.delayed(stl_block)(lats_block, lons_block),
                                                shape=(3, len(dates_weekly), lats_block.size, lons_block.size),
                                                dtype=np.float32)
                blocks.append(block)
                del block
            blocks_total.append(blocks)
            del blocks
        models = dask.array.block(blocks_total)
        del blocks_total

        # transform to separate variables
        coords = {'date': dt_weekly.values, dim1: data[dim1], dim2: data[dim2]}
        # transform to separate variables variables returned from Stack.stl() function
        varnames = ['trend', 'seasonal', 'resid']
        keys_vars = {varname: xr.DataArray(models[varidx], coords=coords) for (varidx, varname) in enumerate(varnames)}
        model = xr.Dataset({**keys_vars})
        del models

        if interactive:
            return model

        self.save_cube(model, name='stl', caption='Seasonal-Trend decomposition using LOESS', chunksize=chunksize)
