#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_dem import SBAS_dem
from .PRM import PRM

class SBAS_stack(SBAS_dem):

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

    # replacement for gmt grdfilter ../topo/dem.grd -D2 -Fg2 -I12s -Gflt.grd
    # use median decimation instead of average
    def get_topo_llt(self, subswath, degrees, geoloc=True, debug=False):
        import xarray as xr
        import numpy as np

        # add buffer around the cropped area for borders interpolation
        dem_area = self.get_dem(subswath, geoloc=geoloc)
        ny = int(np.round(degrees/dem_area.lat.diff('lat')[0]))
        nx = int(np.round(degrees/dem_area.lon.diff('lon')[0]))
        if debug:
            print ('DEBUG: DEM decimation','ny', ny, 'nx', nx)
        dem_area = dem_area.coarsen({'lat': ny, 'lon': nx}, boundary='pad').mean()

        lats, lons, z = xr.broadcast(dem_area.lat, dem_area.lon, dem_area)
        topo_llt = np.column_stack([lons.values.ravel(), lats.values.ravel(), z.values.ravel()])
        return topo_llt

    # aligning for master image
    def stack_ref(self, subswath, debug=False):
        import xarray as xr
        import numpy as np
        import os

    #        err, warn = self.validate()
    #        #print ('err, warn', err, warn)
    #        assert not err and not warn, 'ERROR: Please fix all the issues listed above to continue'

        master_line = list(self.get_master(subswath).itertuples())[0]
        #print (master_line)

        # for master image
        multistem, stem = self.multistem_stem(subswath, master_line.datetime)
        path_stem = os.path.join(self.basedir, stem)
        path_multistem = os.path.join(self.basedir, multistem)

        # generate PRM, LED, SLC
        self.make_s1a_tops(subswath, debug=debug)

        PRM.from_file(path_stem + '.PRM')\
            .set(input_file = path_multistem + '.raw')\
            .update(path_multistem + '.PRM', safe=True)

        self.ext_orb_s1a(subswath, multistem, debug=debug)

        # recalculate after ext_orb_s1a
        earth_radius = PRM.from_file(path_multistem + '.PRM')\
            .calc_dop_orb(inplace=True).update().get('earth_radius')

    # aligning for secondary image
    def stack_rep(self, subswath, date=None, degrees=12.0/3600, debug=False):
        import xarray as xr
        import numpy as np
        import os
        
        # temporary filenames to be removed
        cleanup = []

        master_line = list(self.get_master(subswath).itertuples())[0]
        multistem, stem = self.multistem_stem(subswath, master_line.datetime)
        #print (master_line)

        # define master image parameters
        master = self.PRM(subswath).sel('earth_radius').set(stem=stem, multistem=multistem)

        # prepare coarse DEM for alignment
        # 12 arc seconds resolution is enough, for SRTM 90m decimation is 4x4
        topo_llt = self.get_topo_llt(subswath, degrees=degrees)
        #topo_llt.shape

        line = list(self.get_aligned(subswath, date).itertuples())[0]
        multistem, stem = self.multistem_stem(subswath, line.datetime)
        #print (line)

        # define relative filenames for PRM
        stem_prm    = os.path.join(self.basedir, stem + '.PRM')
        mstem_prm   = os.path.join(self.basedir, multistem + '.PRM')
        master_prm  = os.path.join(self.basedir, master.get("stem") + '.PRM')
        mmaster_prm = os.path.join(self.basedir, master.get("multistem") + '.PRM')

        # TODO: define 1st image for line, in the example we have no more
        tmp_da = 0

        # generate PRM, LED
        self.make_s1a_tops(subswath, date, debug=debug)

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
        #if tmp_da == 0:
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
        if abs(tmp_da) >= 1000:
            prf = tmp_prm.get('PRF')
            tmp_prm.set(tmp_prm.sel('clock_start' ,'clock_stop', 'SC_clock_start', 'SC_clock_stop') - tmp_da/prf/86400.0)
            #raise Exception('TODO: Modifying master PRM by $tmp_da lines...')

        # tmp.PRM defined above from {master}.PRM
        prm1 = tmp_prm.calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug)
        tmpm_dat = prm1.SAT_llt2rat(coords=topo_llt, precise=1, debug=debug)
        prm2 = PRM.from_file(stem_prm).calc_dop_orb(master.get('earth_radius'), inplace=True, debug=debug)
        tmp1_dat = prm2.SAT_llt2rat(coords=topo_llt, precise=1, debug=debug)

        # get r, dr, a, da, SNR table to be used by fitoffset.csh
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
        r_grd_filename = stem_prm[:-4]+'_r.grd'
        r_grd.to_netcdf(r_grd_filename, engine=self.engine)
        # drop the temporary file at the end of the function
        cleanup.append(r_grd_filename)

        a_grd = self.offset2shift(a_xyz, rmax, amax)
        a_grd_filename = stem_prm[:-4]+'_a.grd'
        a_grd.to_netcdf(a_grd_filename, engine=self.engine)
        # drop the temporary file at the end of the function
        cleanup.append(a_grd_filename)

        # generate the image with point-by-point shifts
        # note: it removes calc_dop_orb parameters from PRM file
        # generate PRM, LED
        self.make_s1a_tops(subswath,
                           date=line.Index, mode=1,
                           rshift_fromfile=f'{stem}_r.grd',
                           ashift_fromfile=f'{stem}_a.grd',
                           debug=debug)

        # need to update shift parameter so stitch_tops will know how to stitch
        PRM.from_file(stem_prm).set(PRM.fitoffset(3, 3, offset_dat)).update()

        # echo stitch images together and get the precise orbit
        # use stitch_tops tmp.stitchlist $stem to merge images

        # the raw file does not exist but it works
        PRM.from_file(stem_prm)\
            .set(input_file = f'{multistem}.raw')\
            .update(mstem_prm, safe=True)

        self.ext_orb_s1a(subswath, multistem, date=line.Index, debug=debug)

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
        
        # cleanup
        for filename in cleanup:
            #if os.path.exists(filename):
            os.remove(filename)

    def stack_parallel(self, dates=None, n_jobs=-1, **kwargs):
        from tqdm.auto import tqdm
        import joblib

        if dates is None:
            dates = list(self.get_aligned().index.unique())

        subswaths = self.get_subswaths()

        # prepare master image
        #self.stack_ref()
        with self.tqdm_joblib(tqdm(desc='Reference', total=len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.stack_ref)(subswath, **kwargs) for subswath in subswaths)

        # prepare secondary images
        with self.tqdm_joblib(tqdm(desc='Aligning', total=len(dates)*len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.stack_rep)(subswath, date, **kwargs) \
                                           for date in dates for subswath in subswaths)
