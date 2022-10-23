#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_base import SBAS_base

class SBAS_orbits(SBAS_base):
    # for precision orbit there is only single orbit per day
    # for approximate orbit 2 and maybe more orbits per day are possible
    # so check orbit file for for each subswath
    def download_orbits(self):
        from eof.download import download_eofs

        df = self.df[self.df['orbitpath'].isna()]
    
        # download all the misssed orbit files
        for record in df.itertuples():
            #print (date, mission)
            orbitpath = download_eofs([record.datetime], [record.mission], save_dir=self.basedir)[0]
            #print ('orbitpath', orbitpath)
            self.df.loc[self.df.datetime == record.datetime,'orbitpath'] = orbitpath
