# ----------------------------------------------------------------------------
# InSAR.dev
# 
# This file is part of the InSAR.dev project: https://InSAR.dev
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from insardev_toolkit import tqdm_joblib, tqdm_dask
from .dataset import dataset

class Stack_base(tqdm_joblib, dataset):

    def get_pairs(self, pairs, dates=False):
        """
        Get pairs as DataFrame and optionally dates array.

        Parameters
        ----------
        pairs : np.ndarray, optional
            An array of pairs. If None, all pairs are considered. Default is None.
        dates : bool, optional
            Whether to return dates array. Default is False.
        name : str, optional
            The name of the phase filter. Default is 'phasefilt'.

        Returns
        -------
        pd.DataFrame or tuple
            A DataFrame of pairs. If dates is True, also returns an array of dates.
        """
        import xarray as xr
        import pandas as pd
        import numpy as np
        from glob import glob

        if isinstance(pairs, pd.DataFrame):
            # workaround for baseline_pairs() output
            pairs = pairs.rename(columns={'ref_date': 'ref', 'rep_date': 'rep'})
        elif isinstance(pairs, (xr.DataArray, xr.Dataset)):
            # pairs = pd.DataFrame({
#                 'ref': pairs.coords['ref'].values,
#                 'rep': pairs.coords['rep'].values
#             })
            refs = pairs.coords['ref'].values
            reps = pairs.coords['rep'].values
            pairs = pd.DataFrame({
                'ref': refs if isinstance(refs, np.ndarray) else [refs],
                'rep': reps if isinstance(reps, np.ndarray) else [reps]
            })
        else:
            # Convert numpy array to DataFrame
            # in case of 1d array with 2 items convert to a single pair
            pairs_2d = [pairs] if np.asarray(pairs).shape == (2,) else pairs
            pairs = pd.DataFrame(pairs_2d, columns=['ref', 'rep'])

        # Convert ref and rep columns to datetime format
        pairs['ref'] = pd.to_datetime(pairs['ref'])
        pairs['rep'] = pd.to_datetime(pairs['rep'])
        pairs['pair'] = [f'{ref} {rep}' for ref, rep in zip(pairs['ref'].dt.date, pairs['rep'].dt.date)]
        # Calculate the duration in days and add it as a new column
        #pairs['duration'] = (pairs['rep'] - pairs['ref']).dt.days

        if dates:
            # pairs is DataFrame
            dates = np.unique(pairs[['ref', 'rep']].astype(str).values.flatten())
            return (pairs, dates)
        return pairs

    def get_pairs_matrix(self, pairs):
        """
        Create a matrix based on interferogram dates and pairs.

        Parameters
        ----------
        pairs : pandas.DataFrame or xarray.DataArray or xarray.Dataset
            DataFrame or DataArray containing interferogram date pairs.
        
        Returns
        -------
        numpy.ndarray
            A matrix with one row for every interferogram and one column for every date.
            Each element in the matrix is a float, with 1 indicating the start date,
            -1 indicating the end date, 0 if the date is covered by the corresponding 
            interferogram timeline, and NaN otherwise.

        """
        import numpy as np
        import pandas as pd

        # also define image capture dates from interferogram date pairs
        pairs, dates = self.get_pairs(pairs, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values

        # here are one row for every interferogram and one column for every date
        matrix = []
        for pair in pairs:
            #mrow = [date>=pair[0] and date<=pair[1] for date in dates]
            mrow = [(-1 if date==pair[0] else (1 if date==pair[1] else (0 if date>pair[0] and date<pair[1] else np.nan))) for date in dates]
            matrix.append(mrow)
        matrix = np.stack(matrix).astype(np.float32)
        return matrix

    @staticmethod
    def phase_to_positive_range(phase):
        """
        Convert phase from the range [-pi, pi] to [0, 2pi].
    
        Parameters
        ----------
        phase : array_like
            Input phase values in the range [-pi, pi].
    
        Returns
        -------
        ndarray
            Phase values converted to the range [0, 2pi].
        
        Examples
        --------
        >>> phase_to_positive_range(np.array([-np.pi, -np.pi/2, np.pi, 2*-np.pi-1e-6, 2*-np.pi]))
        array([3.14159265, 4.71238898, 3.14159265, 6.28318431, 0.        ])
        """
        import numpy as np
        return (phase + 2 * np.pi) % (2 * np.pi)
    
    @staticmethod
    def phase_to_symmetric_range(phase):
        """
        Convert phase from the range [0, 2pi] to [-pi, pi].
    
        Parameters
        ----------
        phase : array_like
            Input phase values in the range [0, 2pi].
    
        Returns
        -------
        ndarray
            Phase values converted to the range [-pi, pi].
        
        Examples
        --------
        >>> phase_to_symmetric_range(np.array([0, np.pi, 3*np.pi/2, 2*np.pi]))
        array([ 0.        ,  3.14159265, -1.57079633,  0.        ])
        """
        import numpy as np
        return (phase + np.pi) % (2 * np.pi) - np.pi
