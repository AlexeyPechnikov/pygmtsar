#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to scan *xml files and orbits and make data.in file like to prep_data_linux.csh & prep_data.csh tools
#import pytest

class S1A_Filelist:

    #def _to_io(self, output=None):
    #    return self.df.reset_index().astype(str).apply(lambda row: (' = ').join(row), axis=1)\
    #        .to_csv(output, header=None, index=None)
    #
    #def to_str(self):
    #    return self._to_io()
    #
    #def __str__(self):
    #    return self.to_str()

    def __repr__(self):
        return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    def __init__(self, dirname, master=None):
        import os
        from glob import glob
        import pandas as pd
        from datetime import datetime
        from dateutil.relativedelta import relativedelta
        oneday = relativedelta(days=1)

        # set master image
        self.master = master

        #print ('__init__')
        orbits = glob(os.path.join(dirname, 'S1A_*.EOF'), recursive=True)
        orbits = pd.DataFrame(orbits, columns=['orbitpath'])
        orbits['orbitfile'] = [os.path.split(file)[-1] for file in orbits['orbitpath']]
        orbits['orbitname'] = [os.path.splitext(name)[0] for name in orbits['orbitfile']]
        orbits['date1'] = [name.split('_')[-2][1:9] for name in orbits['orbitname']]
        orbits['date2'] = [name.split('_')[-1][:8] for name in orbits['orbitname']]
        #print (orbits)
    
        metas = glob(os.path.join(dirname, 's1a-iw*.xml'), recursive=True)
        metas = pd.DataFrame(metas, columns=['metapath'])
        metas['metafile'] = [os.path.split(file)[-1] for file in metas['metapath']]
        metas['filename'] = [os.path.splitext(file)[0] for file in metas['metafile']]
        dates = [name[15:30] for name in metas['metafile']]
        dates = [datetime.strptime(date, "%Y%m%dt%H%M%S") for date in dates]
        #print (dates)
        metas['date'] = [date.strftime("%Y-%m-%d") for date in dates]
        metas['date1'] = [(date-oneday).strftime("%Y%m%d") for date in dates]
        metas['date2'] = [(date+oneday).strftime("%Y%m%d") for date in dates]
        #print (filenames)
        #metanames['data'] = [filename[:-4]+'.tiff' for metaname in metanames['metaname']]
        #metanames['dataname'] = [basename[:-4]+'.tiff' for basename in metanames['basename']]
        # TODO: replace F1 to iw*
        metas['stem'] = [f'S1_{date.strftime("%Y%m%d_%H%M%S")}_F1' for date in dates]
        metas['multistem'] = [f'S1_{date.strftime("%Y%m%d")}_ALL_F1' for date in dates]

        datas = glob(os.path.join(dirname, 's1a-iw*.tiff'), recursive=True)
        datas = pd.DataFrame(datas, columns=['datapath'])
        datas['datafile'] = [os.path.split(file)[-1] for file in datas['datapath']]
        datas['filename'] = [os.path.splitext(file)[0] for file in datas['datafile']]
        #print (datas)
    
        self.df = pd.merge(metas, orbits,  how='left', left_on=['date1','date2'], right_on = ['date1','date2'])
        self.df = pd.merge(self.df, datas,  how='left', left_on=['filename'], right_on = ['filename']).set_index('date')
        del self.df['date1']
        del self.df['date2']

    def get_master(self):
        if self.master is None:
            raise Exception('Set master image first')
        idx = self.master
        return self.df.loc[[idx]]

    def get_aligned(self):
        if self.master is None:
            raise Exception('Set master image first')
        idx = self.df.index.difference([self.master])
        return self.df.loc[idx]

    def to_file(self, filename):
        """
        Save data.in file like to prep_data_linux.csh & prep_data.csh tools
        """
        if self.master is None:
            raise Exception('Set master image first')
        line = '\n'.join(self.get_master().apply(lambda row: f'{row.filename}:{row.orbitfile}', axis=1).values)
        lines = '\n'.join(self.get_aligned().apply(lambda row: f'{row.filename}:{row.orbitfile}', axis=1).values)
        with open(filename, 'wt') as f:
            f.write(line+'\n'+lines+'\n')
        return self

    def set_master(self, master):
        self.master = master
        return self

    def init(self, output_directory, topo_filename):
        """
        Initialize directory for processing
        """
        import os
        import shutil

        shutil.rmtree(output_directory, ignore_errors=True)
        os.makedirs(output_directory, exist_ok=True)

        metas     = list(self.df['metapath'].values)
        datas = list(self.df['datapath'].values)
        orbits    = list(self.df['orbitpath'].values)
        for file in metas + datas + orbits:
            path = os.path.join(output_directory, os.path.split(file)[-1])
            print (f'{file} -> {path}')
            os.symlink(os.path.relpath(file,output_directory), path)

        path = os.path.join(output_directory, 'dem.grd')
        os.symlink(os.path.relpath(topo_filename,output_directory), path)

        return self




#filelist = S1A_Filelist('raw_orig').set_master(MASTER)
#filelist.df
#filelist.get_master()
#filelist.get_aligned()
#filelist.to_file('raw/data.in.new')
