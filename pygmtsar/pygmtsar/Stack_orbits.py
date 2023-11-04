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
    def download_orbits(self, strict=False):
        """
        Download missed orbits for all the Stack scenes.

        Returns
        -------
        None

        Examples
        --------
        stack.download_orbits()
        """
        #from eof.download import download_eofs

        print ('NOTICE: this function is removed due to continuous downloading issues from free access servers.')
        print ('NOTICE: use ASF module to download scenes and orbits.')
        print ('NOTICE: or call "eof" command-line tool in your data directory.')

#         if not strict:
#             # not strict (fast) mode — download a single orbit file for (mission, date) pairs
#             # for all (date, mission) unique pairs download the corresponding orbit files
#             df = self.df[self.df['orbitpath'].isna()][['datetime','mission']].reset_index()\
#                 .drop_duplicates(['date','mission']).sort_values('datetime')
#             if len(df) == 0:
#                 return
#             orbitpaths = download_eofs(df.datetime.tolist(), df.mission.tolist(), save_dir=self.basedir)
#             # dataframe is already sorted
#             df['orbitpath'] = sorted(orbitpaths, key=lambda x: x[-19:])
#             for record in df.itertuples():
#                 self.df.loc[(self.df.index == record.date)&(self.df.mission == record.mission),'orbitpath'] = record.orbitpath
#         else:
#             # strict mode — iterate all the records and download all the orbits separately
#             df = self.df[self.df['orbitpath'].isna()]
#             # download all the misssed orbit files
#             for record in df.itertuples():
#                 #print (date, mission)
#                 orbitpath = download_eofs([record.datetime], [record.mission], save_dir=self.basedir)[0]
#                 #print ('orbitpath', orbitpath)
#                 self.df.loc[self.df.datetime == record.datetime,'orbitpath'] = orbitpath
