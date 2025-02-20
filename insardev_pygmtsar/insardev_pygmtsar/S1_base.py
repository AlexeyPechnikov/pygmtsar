# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from insardev_toolkit import tqdm_joblib, tqdm_dask
from .dataset import dataset

class S1_base(tqdm_joblib, dataset):

    def __repr__(self):
        return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    def to_dataframe(self):
        """
        Return a Pandas DataFrame for all Stack scenes.

        Returns
        -------
        pandas.DataFrame
            The DataFrame containing Stack scenes.

        Examples
        --------
        df = stack.to_dataframe()
        """
        return self.df

    def get_prefix(self, burst, date=None):
        import os

        assert date is None or date is False or len(date)==10, f'ERROR: date format is not yyyy-mm-dd (burst={burst} date={date})'
        # TODO
        assert len(burst)!=10, 'ERROR: mixed burst and date arguments (burst={burst} date={date})'

        path = os.path.join(self.basedir, burst)
        if not os.path.isdir(path):
            os.makedirs(path)

        if date == False:
            return os.path.join(burst, '')
            
        # use reference datetime if not defined
#         if date is None or date  == self.reference:
#             df = self.get_reference(burst)
#         else:
#             df = self.get_repeat(burst, date)    
#         name = df.burst.iloc[0]
        if date is None:
            date = self.reference
        name = date.replace('-','')
        return os.path.join(burst, name)

    def set_reference(self, reference):
        """
        Define reference scene for Stack object.

        Parameters
        ----------
        reference : str
            Date string representing the reference scene.

        Returns
        -------
        Stack
            Modified instance of the Stack class.

        Examples
        --------
        Set the reference scene to '2022-01-20':
        stack.set_reference('2022-01-20')
        """
        if reference is None:
            #print ('NOTE: reference scene is None, Stack.set_reference() command is ignored')
            if self.reference is None:
                self.reference = self.df.index.get_level_values(1)[0]
                print (f'NOTE: auto set reference scene {self.reference}. You can change it like Stack.set_reference("{self.reference}")')
            return self
        assert reference in self.df.index.get_level_values(1), f'Reference scene not found: {reference}'
        self.reference = reference
        return self

    def get_reference(self, burst):
        """
        Return dataframe reference record.

        Parameters
        ----------
        None

        Returns
        -------
        pd.DataFrame
            The DataFrame containing reference record.
        """
        df = self.df.loc[[(burst, self.reference)]]
        assert len(df) > 0, f'Reference record not found'
        return df
        
    def get_repeat(self, burst, date=None):
        """
        Return dataframe repeat records (excluding reference).

        Parameters
        ----------
        date : datetime, optional
            The date for which to return repeat records. If None, all dates are considered. Default is None.

        Returns
        -------
        pd.DataFrame
            The DataFrame containing repeat records for the specified date.
        """
        if date is None:
            df_filtered = self.df[self.df.index.get_level_values(0) == burst]
            idx_reference = self.df.index[self.df.index.get_level_values(1) == self.reference]
            return df_filtered.loc[df_filtered.index.difference(idx_reference)]

        assert not date == self.reference, f'ERROR: repeat date cannot be equal to reference date "{date}"'
        return self.df.loc[[(burst, date)]]
