# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_prm import Stack_prm

class Stack_orbits(Stack_prm):
    # for precision orbit there is only single orbit per day
    # for approximate orbit 2 and maybe more orbits per day are possible
    def download_orbits(self):
        """
        Download missed orbits for all the Stack scenes.

        Returns
        -------
        None


        Examples
        --------
        stack.download_orbits()
        """
        from eof.download import download_eofs

        # for all (date, mission) unique pairs download the corresponding orbit files
        df = self.df[self.df['orbitpath'].isna()][['mission']].reset_index().drop_duplicates().sort_values('date')
        if len(df) == 0:
            return
        orbitpaths = download_eofs(df.date.tolist(), df.mission.tolist(), save_dir=self.basedir)
        # dataframe is already sorted
        df['orbitpath'] = sorted(orbitpaths, key=lambda x: x[-19:])
        for record in df.itertuples():
            self.df.loc[(self.df.index == record.date)&(self.df.mission == record.mission),'orbitpath'] = record.orbitpath
