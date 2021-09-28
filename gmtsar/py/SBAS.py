#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to scan *xml files and orbits and make data.in file like to prep_data_linux.csh & prep_data.csh tools
#import pytest

class SBAS:
    import contextlib
    @staticmethod
    @contextlib.contextmanager
    def tqdm_joblib(tqdm_object):
        import joblib
        """Context manager to patch joblib to report into tqdm progress bar given as argument"""
        class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)

            def __call__(self, *args, **kwargs):
                tqdm_object.update(n=self.batch_size)
                return super().__call__(*args, **kwargs)

        old_batch_callback = joblib.parallel.BatchCompletionCallBack
        joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
        try:
            yield tqdm_object
        finally:
            joblib.parallel.BatchCompletionCallBack = old_batch_callback
            tqdm_object.close()

    @staticmethod
    def nearest_grid(in_grid,  search_radius_pixels=300):
        from scipy.spatial import cKDTree
        import xarray as xr
        import numpy as np

        ys, xs = np.meshgrid(range(in_grid.shape[1]), range(in_grid.shape[0]))
        ys = ys.reshape(-1)
        xs = xs.reshape(-1)
        zs = in_grid.values.reshape(-1)
        mask = np.where(~np.isnan(zs))
        # on regular source grid some options should be redefined for better performance
        tree = cKDTree(np.column_stack((ys[mask],xs[mask])), compact_nodes=False, balanced_tree=False)
        # use distance_limit
        d, inds = tree.query(np.column_stack((ys,xs)), k = 1, distance_upper_bound=search_radius_pixels, workers=8)
        # replace not available indexes by zero (see distance_upper_bound)
        fakeinds = np.where(~np.isinf(d), inds, 0)
        # produce the same output array as dataset to be able to add global attributes
        values = np.where(~np.isinf(d), zs[mask][fakeinds], np.nan).reshape(in_grid.shape)
        out_grid = xr.DataArray(values, coords=in_grid.coords, name='z')
        return out_grid

    def __repr__(self):
        return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    def __init__(self, datadir, dem_filename, basedir):
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
        
        # processing directory
        self.basedir = basedir

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

    def get_aligned(self, date=None):
        if self.master is None:
            raise Exception('Set master image first')
        #if self.master is not None and self.master == date:
        #    raise Exception('Requested image is master image')

        if date is None:
            idx = self.df.index.difference([self.master])
        else:
            idx = [date]
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
    
    # buffer required to get correct (binary) results from SAT_llt2rat tool
    # minimum buffer size: 8 arc seconds for 90 m DEM
    def get_dem(self, geoloc=False, buffer_degrees = 12./3600):
        import xarray as xr
                
        dem = xr.open_dataarray(self.dem_filename)
        if geoloc is False:
            return dem
        
        geoloc = self.geoloc()
        ymin, ymax = geoloc.latitude.min(), geoloc.latitude.max()
        #print ('ymin, ymax', ymin, ymax)
        xmin, xmax = geoloc.longitude.min(), geoloc.longitude.max()
        #print ('xmin, xmax', xmin, xmax)
        return dem.sel(lat=slice(ymin-buffer_degrees, ymax+buffer_degrees),
                       lon=slice(xmin-buffer_degrees, xmax+buffer_degrees))

    def topo_ra(self, method='cubic'):
        import xarray as xr
        import os

        trans_dat_file = os.path.join(self.basedir, 'trans.dat')
        topo_ra_file = os.path.join(self.basedir, 'topo_ra.grd')

        if os.path.exists(topo_ra_file):
            os.remove(topo_ra_file)

        self.PRM().topo_ra(self.get_dem(geoloc=True),
                           trans_dat_tofile=trans_dat_file,
                           topo_ra_tofile=topo_ra_file,
                           method=method)

    # replacement for gmt grdfilter ../topo/dem.grd -D2 -Fg2 -I12s -Gflt.grd
    def get_topo_llt(self, degrees, geoloc=True):
        import xarray as xr
        import numpy as np

        # add buffer around the cropped area for borders interpolation
        dem_area = self.get_dem(geoloc=geoloc, buffer_degrees=2*degrees)
        ny = int(np.round(degrees/dem_area.lat.diff('lat')[0]))
        nx = int(np.round(degrees/dem_area.lon.diff('lon')[0]))
        #print ('DEM decimation','ny', ny, 'nx', nx)
        dem_area = dem_area.coarsen({'lat': ny, 'lon': nx}, boundary='pad').mean()

        lats, lons, z = xr.broadcast(dem_area.lat, dem_area.lon, dem_area)
        topo_llt = np.column_stack([lons.values.ravel(), lats.values.ravel(), z.values.ravel()])
        return topo_llt

    def offset2shift(self, xyz, rmax, amax, method='linear'):
        import xarray as xr
        import numpy as np
        from scipy.interpolate import griddata

        # use center pixel GMT registration mode
        rngs = np.arange(8/2, rmax+8/2, 8)
        azis = np.arange(4/2, amax+4/2, 4)
        grid_r, grid_a = np.meshgrid(rngs, azis)

        grid = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (grid_r, grid_a), method=method)
        da = xr.DataArray(np.flipud(grid), coords={'y': azis, 'x': rngs}, name='z')
        return da

    def to_dataframe(self):
        return self.df

    def ext_orb_s1a(self, stem, date=None, debug=False):
        import os
        import subprocess

        if date is None or date == self.master:
            df = self.get_master()
        else:
            df = self.df.loc[[date]]

        orbit = os.path.relpath(df['orbitpath'][0], self.basedir)
    
        argv = ['ext_orb_s1a', f'{stem}.PRM', orbit, stem]
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('ext_orb_s1a', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('ext_orb_s1a', stdout_data.decode('ascii'))

        return
    
    # produce LED and PRM in basedir
    # when date=None work on master image
    def make_s1a_tops(self, date=None, mode=0, rshift_fromfile=None, ashift_fromfile=None, debug=False):
        import os
        import subprocess

        if date is None or date == self.master:
            date = self.master
            # for master image mode should be 1
            mode = 1
        df = self.df.loc[[date]]
    
        xmlfile = os.path.relpath(df['metapath'][0], self.basedir)
        datafile = os.path.relpath(df['datapath'][0], self.basedir)
        #orbit = os.path.relpath(df['orbitfile'][0], self.basedir)
        stem = df['stem'][0]
    
        argv = ['make_s1a_tops', xmlfile, datafile, stem, str(mode)]
        if rshift_fromfile is not None:
            argv.append(rshift_fromfile)
        if ashift_fromfile is not None:
            argv.append(ashift_fromfile)
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('make_s1a_tops', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('make_s1a_tops', stdout_data.decode('ascii'))
    
        self.ext_orb_s1a(stem, date, debug=debug)
    
        return

    # aligning for master image
    def stack_ref(self, debug=False):
        import xarray as xr
        import numpy as np
        import os
        from PRM import PRM

        if self.master is None:
            raise Exception('Set master image first')
                
        master_line = list(self.get_master().itertuples())[0]
        #print (master_line)

        # for master image
        path_stem = os.path.join(self.basedir, master_line.stem)
        path_multistem = os.path.join(self.basedir, master_line.multistem)

        # generate prms and leds
        self.make_s1a_tops(debug=debug)

        PRM.from_file(path_stem + '.PRM')\
            .set(input_file = path_multistem + '.raw')\
            .update(path_multistem + '.PRM', safe=True)

        self.ext_orb_s1a(master_line.multistem, debug=debug)

        # recalculate after ext_orb_s1a
        earth_radius = PRM.from_file(path_multistem + '.PRM')\
            .calc_dop_orb(inplace=True).update().get('earth_radius')

    # aligning for secondary image
    def stack_rep(self, date=None, degrees=12.0/3600, debug=False):
        import xarray as xr
        import numpy as np
        import os
        from PRM import PRM

        if self.master is None:
            raise Exception('Set master image first')

        master_line = list(self.get_master().itertuples())[0]
        #print (master_line)

        # define master image parameters
        master = self.PRM().sel('earth_radius').set(stem=master_line.stem, multistem=master_line.multistem)

        # prepare coarse DEM for alignment
        # 12 arc seconds resolution is enough, for SRTM 90m decimation is 4x4
        topo_llt = self.get_topo_llt(degrees=degrees)
        #topo_llt.shape

        line = list(self.get_aligned(date).itertuples())[0]
        #print (line)

        # define relative filenames for PRM
        stem_prm    = os.path.join(self.basedir, line.stem + '.PRM')
        mstem_prm   = os.path.join(self.basedir, line.multistem + '.PRM')
        master_prm  = os.path.join(self.basedir, master.get("stem") + '.PRM')
        mmaster_prm = os.path.join(self.basedir, master.get("multistem") + '.PRM')

        # TODO: define 1st image for line, in the example we have no more
        tmp_da = 0

        # generate prms and leds
        self.make_s1a_tops(date, debug=debug)

        # compute the time difference between first frame and the rest frames
        t1, prf = PRM.from_file(stem_prm).get('clock_start', 'PRF')
        t2      = PRM.from_file(stem_prm).get('clock_start')
        nl = int((t2 - t1)*prf*86400.0+0.2)
        #echo "Shifting the master PRM by $nl lines..."

        # Shifting the master PRM by $nl lines...
        # shift the super-masters PRM based on $nl so SAT_llt2rat gives precise estimate
        prm1 = PRM.from_file(master_prm)
        prm1.set(prm1.sel('clock_start' ,'clock_stop', 'SC_clock_start', 'SC_clock_stop') + nl/prf/86400.0)
        tmp_prm = prm1

        # compute whether there are any image offset
        if tmp_da == 0:
            # tmp_prm defined above from {master}.PRM
            prm1 = tmp_prm.calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug)
            prm2 = PRM.from_file(stem_prm).calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug).update()
            lontie,lattie = prm1.SAT_baseline(prm2, debug=debug).get('lon_tie_point', 'lat_tie_point')
            tmp_am = prm1.SAT_llt2rat(coords=[lontie, lattie, 0], precise=1, debug=debug)[1]
            tmp_as = prm2.SAT_llt2rat(coords=[lontie, lattie, 0], precise=1, debug=debug)[1]
            # bursts look equal to rounded result int(np.round(...))
            tmp_da = int(tmp_as - tmp_am)
            #print ('tmp_am', tmp_am, 'tmp_as', tmp_as, 'tmp_da', tmp_da)

        # in case the images are offset by more than a burst, shift the super-master's PRM again
        # so SAT_llt2rat gives precise estimate
        if abs(tmp_da) < 1000:
            # tmp.PRM defined above from {master}.PRM
            prm1 = tmp_prm.calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug)
            tmpm_dat = prm1.SAT_llt2rat(coords=topo_llt, precise=1, debug=debug)
            prm2 = PRM.from_file(stem_prm).calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug)
            tmp1_dat = prm2.SAT_llt2rat(coords=topo_llt, precise=1, debug=debug)
        else:
            raise Exception('TODO: Modifying master PRM by $tmp_da lines...')

        # echo get r, dr, a, da, SNR table to be used by fitoffset.csh
        offset_dat0 = np.hstack([tmpm_dat, tmp1_dat])
        func = lambda row: [row[0],row[5]-row[0],row[1],row[6]-row[1],100]
        offset_dat = np.apply_along_axis(func, 1, offset_dat0)

        # define radar coordinates extent
        rmax, amax = PRM.from_file(stem_prm).get('num_rng_bins','num_lines')

        # prepare the offset parameters for the stitched image
        # set the exact borders in radar coordinates
        par_tmp = offset_dat[(offset_dat[:,0]>0) & (offset_dat[:,0]<rmax) & (offset_dat[:,2]>0) & (offset_dat[:,2]<amax)]
        par_tmp[:,2] += nl
        if abs(tmp_da) >= 1000:
            par_tmp[:,2] -= tmp_da
            par_tmp[:,3] += tmp_da

        # prepare the rshift and ashift look up table to be used by make_s1a_tops
        # use tmp_dat instead of offset_dat
        r_xyz = offset_dat[:,[0,2,1]]
        a_xyz = offset_dat[:,[0,2,3]]

        r_grd = self.offset2shift(r_xyz, rmax, amax)
        r_grd.to_netcdf(stem_prm[:-4]+'_r.grd')

        a_grd = self.offset2shift(a_xyz, rmax, amax)
        a_grd.to_netcdf(stem_prm[:-4]+'_a.grd')

        # generate the image with point-by-point shifts
        # note: it removes calc_dop_orb parameters from PRM file
        self.make_s1a_tops(date=line.Index, mode=1,
                           rshift_fromfile=f'{line.stem}_r.grd',
                           ashift_fromfile=f'{line.stem}_a.grd',
                           debug=debug)

        # need to update shift parameter so stitch_tops will know how to stitch
        PRM.from_file(stem_prm).set(PRM.fitoffset(3, 3, offset_dat)).update()

        # echo stitch images together and get the precise orbit
        # use stitch_tops tmp.stitchlist $stem to merge images

        # the raw file does not exist but it works
        PRM.from_file(stem_prm)\
            .set(input_file = f'{line.multistem}.raw')\
            .update(mstem_prm, safe=True)

        self.ext_orb_s1a(line.multistem, date=line.Index, debug=debug)

        # Restoring $tmp_da lines shift to the image... 
        PRM.from_file(mstem_prm).set(ashift=0 if abs(tmp_da) < 1000 else tmp_da, rshift=0).update()

        # that is safe to rewrite source files
        prm1 = PRM.from_file(mmaster_prm)
        prm1.resamp(PRM.from_file(mstem_prm),
                    alignedSLC_tofile=mstem_prm[:-4]+'.SLC',
                    interp=1, debug=debug
        ).to_file(mstem_prm)

        PRM.from_file(mstem_prm).set(PRM.fitoffset(3, 3, par_tmp)).update()

        PRM.from_file(mstem_prm).calc_dop_orb(master.get('earth_radius'), 0, inplace=True, debug=debug).update()

    def stack_parallel(self, dates=None, **kwargs):
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib
    
        if dates is None:
            if self.master is None:
                raise Exception('Set master image or define argument dates')
            dates = list(self.get_aligned().index)

        # prepare master image
        self.stack_ref()

        # prepare secondary images
        with self.tqdm_joblib(notebook.tqdm(desc="Aligning", total=len(dates))) as progress_bar:
            joblib.Parallel(n_jobs=-1)(joblib.delayed(self.stack_rep)(date, **kwargs) for date in dates)

    def intf(self, pair, geocode=True, **kwargs):
        from PRM import PRM
        import os

        # extract dates from pair
        date1, date2 = pair

        prm_ref = self.PRM(date1)
        prm_rep = self.PRM(date2)

        # TODO: define function for geocoding
        # apply the user function AFTER geocoding
        if geocode:
            func = None
            if 'func' in kwargs:
                func = kwargs['func']
            def func_geocoding(da_ra):
                da_ll = da_ra
                if func is not None:
                    return func(da_ll)
                return da_ll
            # reassign callback function
            kwargs['func'] = func_geocoding

        print ('SBAS intf kwargs', kwargs)
        prm_ref.intf(prm_rep, basedir=self.basedir, **kwargs)

    def intf_parallel(self, pairs, **kwargs):
        import pandas as pd
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        with self.tqdm_joblib(notebook.tqdm(desc="Interferograms", total=len(pairs))) as progress_bar:
            joblib.Parallel(n_jobs=-1)(joblib.delayed(self.intf)(pair, **kwargs) for pair in pairs)

    def baseline_table(self):
        import pandas as pd
        
        prm_ref = self.PRM()
        data = []
        for date in self.df.index:
            prm_rep = self.PRM(date)
            ST0 = prm_rep.get('SC_clock_start')
            DAY = int(ST0 % 1000)
            YR = int(ST0/1000) - 2014
            YDAY = YR * 365 + DAY
            #print (f'YR={YR}, DAY={DAY}')
            BPL, BPR = prm_ref.SAT_baseline(prm_rep).get('B_parallel', 'B_perpendicular')
            data.append({'date':date, 'ST0':ST0, 'YDAY':YDAY, 'BPL':BPL, 'BPR':BPR})
            #print (date, ST0, YDAY, BPL, BPR)
        return pd.DataFrame(data).set_index('date')

    def baseline_pairs(self, days=100, meters=150):   
        import numpy as np
        import pandas as pd
     
        tbl = self.baseline_table()
        data = []
        for line1 in tbl.itertuples():
        #for line1 in tbl.loc[['2015-01-21']].itertuples():
            for line2 in tbl.itertuples():
            #for line2 in tbl.loc[['2015-03-10']].itertuples():
                #print (line1, line2)
                if not (line1.YDAY < line2.YDAY and line2.YDAY - line1.YDAY < days):
                    continue
                if not (abs(line1.BPR - line2.BPR)< meters):
                    continue
            
                data.append({'ref_date':line1.Index, 'rep_date': line2.Index,
                             'ref_timeline': np.round(line1.YDAY/365.25+2014, 2), 'ref_baseline': np.round(line1.BPR, 2),
                             'rep_timeline': np.round(line2.YDAY/365.25+2014, 2), 'rep_baseline': np.round(line2.BPR, 2)})
            
                #print (line1.Index, line2.Index,
                #       np.round(line1.YDAY/365.25+2014, 2), np.round(line1.BPR, 2),
                #       np.round(line2.YDAY/365.25+2014, 2), np.round(line2.BPR, 2))
        return pd.DataFrame(data)

    # TODO
    def baseline_triplets(self, days, meters):
        return

    def PRM(self, date=None, multi=True):
        from PRM import PRM
        import os

        if date is None or date == self.master:
            line = self.get_master()
        else:
            line = self.get_aligned(date)
        #print (line)
        if multi:
            stem = line.multistem[0]
        else:
            stem = line.stem[0]
        filename = os.path.join(self.basedir, f'{stem}.PRM')
        #print (filename)
        return PRM.from_file(filename)

    def unwrap_parallel(self, pairs, **kwargs):
        import pandas as pd
        #from tqdm import tqdm
        from tqdm import notebook
        import joblib

        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        with self.tqdm_joblib(notebook.tqdm(desc="Unwrapping", total=len(pairs))) as progress_bar:
            joblib.Parallel(n_jobs=-1)(joblib.delayed(self.unwrap)(pair, **kwargs) for pair in pairs)


    # -s for SMOOTH mode and -d for DEFO mode when DEFOMAX_CYCLE should be defined in the configuration
    # DEFO mode (-d) and DEFOMAX_CYCLE=0 is equal to SMOOTH mode (-s)
    def unwrap(self, pair, threshold=0.1, conf=None, func=None, debug=False):
        import xarray as xr
        import numpy as np
        import os
        import subprocess

        # save to NetCDF grid
        compression = dict(zlib=True, complevel=3, chunksizes=[128,128])
    
        if conf is None:
            conf = self.PRM().snaphu_config(defomax=0)

        # extract dates from pair
        date1, date2 = pair

        basename = os.path.join(self.basedir, f'{date1}_{date2}_').replace('-','')
        #print ('basename', basename)

        # invalid pixels mask
        mask = xr.open_dataarray(basename + 'mask.grd')
    
        phase = xr.open_dataarray(basename + 'phasefilt.grd')
        phase_in = basename + 'unwrap.phase'
        phase.where(~np.isnan(phase),0).values.tofile(phase_in)

        corr = xr.open_dataarray(basename + 'corr.grd')
        # mask invalid pixels on correlation matrix
        corr = (corr*mask).where(corr>=threshold)
        corr_in = basename + 'unwrap.corr'
        corr.where(~np.isnan(corr),0).values.tofile(corr_in)

        unwrap_out = basename + 'unwrap.out'

        argv = ['snaphu', phase_in, str(phase.shape[1]), '-c', corr_in,
                '-f', '/dev/stdin', '-o', unwrap_out, '-d']
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='ascii', bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=conf)
        if len(stderr_data) > 0 and debug:
            print ('snaphu', stderr_data)

        if debug:
            print ('snaphu', stdout_data)

        # read results
        values = np.fromfile(unwrap_out, dtype=np.float32).reshape(phase.shape)
        #values = np.frombuffer(stdout_data, dtype=np.float32).reshape(mask.shape)
        # mask invalid pixels on unwrapped results
        unwrap = mask * xr.DataArray(values, phase.coords, name='z')
        # apply user-defined function for post-processing
        if func is not None:
            unwrap = func(corr, unwrap)
        if os.path.exists(basename + 'unwrap.grd'):
            os.remove(basename + 'unwrap.grd')
        unwrap.to_netcdf(basename + 'unwrap.grd', encoding={'z': compression})

        for tmp_file in [phase_in, corr_in, unwrap_out]:
            #print ('tmp_file', tmp_file)
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

#filelist = SBAS('raw_orig').set_master(MASTER)
#filelist.df
#filelist.get_master()
#filelist.get_aligned()
#filelist.to_file('raw/data.in.new')
