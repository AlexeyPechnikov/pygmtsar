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
    # so check orbit file for for each subswath
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

        df = self.df[self.df['orbitpath'].isna()]
    
        # download all the misssed orbit files
        for record in df.itertuples():
            #print (date, mission)
            orbitpath = download_eofs([record.datetime], [record.mission], save_dir=self.basedir)[0]
            #print ('orbitpath', orbitpath)
            self.df.loc[self.df.datetime == record.datetime,'orbitpath'] = orbitpath
