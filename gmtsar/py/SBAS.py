#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to scan *xml files and orbits and make data.in file like to prep_data_linux.csh & prep_data.csh tools
#import pytest

class SBAS:

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

    def __init__(self, datadir, dem_filename):
        import os
        from glob import glob
        import pandas as pd
        from datetime import datetime
        from dateutil.relativedelta import relativedelta
        oneday = relativedelta(days=1)

        self.dem_filename = os.path.relpath(dem_filename,'.')
        #print ('dem_filename', self.dem_filename)

        # master image
        self.master = None

        orbits = glob(os.path.join(datadir, 'S1?*.EOF'), recursive=True)
        orbits = pd.DataFrame(orbits, columns=['orbitpath'])
        orbits['orbitfile'] = [os.path.split(file)[-1] for file in orbits['orbitpath']]
        orbits['orbitname'] = [os.path.splitext(name)[0] for name in orbits['orbitfile']]
        orbits['date1'] = [name.split('_')[-2][1:9] for name in orbits['orbitname']]
        orbits['date2'] = [name.split('_')[-1][:8] for name in orbits['orbitname']]
        #print (orbits)

        metas = glob(os.path.join(datadir, 's1?-iw*.xml'), recursive=True)
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

        datas = glob(os.path.join(datadir, 's1?-iw*.tiff'), recursive=True)
        datas = pd.DataFrame(datas, columns=['datapath'])
        datas['datafile'] = [os.path.split(file)[-1] for file in datas['datapath']]
        datas['filename'] = [os.path.splitext(file)[0] for file in datas['datafile']]
        #print (datas)

        self.df = pd.merge(metas, orbits,  how='left', left_on=['date1','date2'], right_on = ['date1','date2'])
        self.df = pd.merge(self.df, datas,  how='left', left_on=['filename'], right_on = ['filename']).set_index('date')
        del self.df['date1']
        del self.df['date2']

    def set_master(self, master):
        self.master = master
        return self

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

    def geoloc(self, date=None):
        import pandas as pd
        import xmltodict

        if date is None:
            if self.master is None:
                raise Exception('Set master image or define argument date')
            date = self.master

        filename = self.df.loc[date,'metapath']
        with open(filename) as fd:
            doc = xmltodict.parse(fd.read())
        #doc['geolocationGrid']
        geoloc = doc['product']['geolocationGrid']['geolocationGridPointList']
        # check data consistency
        assert int(geoloc['@count']) == len(geoloc['geolocationGridPoint'])
        geoloc_df = pd.DataFrame(geoloc['geolocationGridPoint']).applymap(lambda val : pd.to_numeric(val,errors='ignore'))
        return geoloc_df
        
    def get_dem(self, geoloc=False):
        import xarray as xr
        
        dem = xr.open_dataarray(self.dem_filename)
        if geoloc is False:
            return dem
        
        geoloc = self.geoloc()
        ymin, ymax = geoloc.latitude.min(), geoloc.latitude.max()
        #print ('ymin, ymax', ymin, ymax)
        xmin, xmax = geoloc.longitude.min(), geoloc.longitude.max()
        #print ('xmin, xmax', xmin, xmax)
        return dem.sel(lat=slice(ymin, ymax), lon=slice(xmin,xmax))

    # replacement for gmt grdfilter ../topo/dem.grd -D2 -Fg2 -I12s -Gflt.grd
    def get_topo_llt(self, degrees, geoloc=True):
        import xarray as xr
        import numpy as np

        dem_area = self.get_dem(geoloc=geoloc)
        ny = int(np.round(degrees/dem_area.lat.diff('lat')[0]))
        nx = int(np.round(degrees/dem_area.lon.diff('lon')[0]))
        #print ('DEM decimation','ny', ny, 'nx', nx)
        dem_area = dem_area.coarsen({'lat': ny, 'lon': nx}, boundary='pad').median()

        lats, lons, z = xr.broadcast(dem_area.lat, dem_area.lon, dem_area)
        topo_llt = np.column_stack([lons.values.ravel(), lats.values.ravel(), z.values.ravel()])
        return topo_llt

    def to_dataframe(self):
        return self.df

    def ext_orb_s1a(self, basedir, stem, date=None):
        import os
        import subprocess

        if date is None or date == self.master:
            df = self.get_master()
        else:
            df = self.df.loc[[date]]

        orbit = os.path.relpath(df['orbitpath'][0], basedir)
        #stem = df['stem'][0]
    
        argv = ['ext_orb_s1a', f'{stem}.PRM', orbit, stem]
        print ('argv', argv)
        p = subprocess.Popen(argv, stderr=subprocess.PIPE, cwd=basedir)
        stderr_data = p.communicate()[1]
        if len(stderr_data) > 0:
            print (stderr_data.decode('ascii'))
        return
    
    # produce LED and PRM in basedir
    # when date=None work on master image
    def make_s1a_tops(self, basedir, date=None, mode=0, rshift_fromfile=None, ashift_fromfile=None):
        import os
        import subprocess

        if date is None or date == self.master:
            date = self.master
            # for master image mode should be 1
            mode = 1
        df = self.df.loc[[date]]
    
        xmlfile = os.path.relpath(df['metapath'][0], basedir)
        datafile = os.path.relpath(df['datapath'][0], basedir)
        #orbit = os.path.relpath(df['orbitfile'][0], basedir)
        stem = df['stem'][0]
        #mmaster = df['multistem'][0]
        #print (file, orbit, master)
    
        # generate prms and leds
        #!cd raw && make_s1a_tops {file}.xml {file}.tiff {master} 1
        #!cd raw && ext_orb_s1a {master}.PRM {orbit} {master}
    
        argv = ['make_s1a_tops', xmlfile, datafile, stem, str(mode)]
        if rshift_fromfile is not None:
            argv.append(rshift_fromfile)
        if ashift_fromfile is not None:
            argv.append(ashift_fromfile)
        print ('argv', argv)
        p = subprocess.Popen(argv, stderr=subprocess.PIPE, cwd=basedir)
        stderr_data = p.communicate()[1]
        if len(stderr_data) > 0:
            print (stderr_data.decode('ascii'))
    
        self.ext_orb_s1a(basedir, stem, date)
    
        return

    def preproc(self, basedir, date=None):
        import os
        from PRM import PRM

        # for master image
        if date is None or date == self.master:
            master_line = list(self.get_master().itertuples())[0]
            #print (master_line)

            path_stem = os.path.join(basedir, master_line.stem)
            path_multistem = os.path.join(basedir, master_line.multistem)

            # generate prms and leds
            self.make_s1a_tops(basedir)

            PRM.from_file(path_stem + '.PRM')\
                .set(input_file = path_multistem + '.raw')\
                .update(path_multistem + '.PRM', safe=True)

            self.ext_orb_s1a(basedir, master_line.multistem)

            # recalculate after ext_orb_s1a
            earth_radius = PRM.from_file(path_multistem + '.PRM')\
                .calc_dop_orb(inplace=True).update().get('earth_radius')

            return PRM().set(earth_radius=earth_radius, stem=master_line.stem, multistem=master_line.multistem)

        # TODO: aligning for secondary images

    def baseline_table(self, days, meters):
        return


#filelist = SBAS('raw_orig').set_master(MASTER)
#filelist.df
#filelist.get_master()
#filelist.get_aligned()
#filelist.to_file('raw/data.in.new')
