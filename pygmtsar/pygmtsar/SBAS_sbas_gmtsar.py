#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_detrend import SBAS_detrend

class SBAS_sbas_gmtsar(SBAS_detrend):

    def sbas(self, baseline_pairs, smooth=0, atm=0, debug=False):
        """
         USAGE: sbas intf.tab scene.tab N S xdim ydim [-atm ni] [-smooth sf] [-wavelength wl] [-incidence inc] [-range -rng] [-rms] [-dem]

         input:
          intf.tab             --  list of unwrapped (filtered) interferograms:
           format:   unwrap.grd  corr.grd  ref_id  rep_id  B_perp
          scene.tab            --  list of the SAR scenes in chronological order
           format:   scene_id   number_of_days
           note:     the number_of_days is relative to a reference date
          N                    --  number of the interferograms
          S                    --  number of the SAR scenes
          xdim and ydim        --  dimension of the interferograms
          -smooth sf           --  smoothing factors, default=0
          -atm ni              --  number of iterations for atmospheric correction, default=0(skip atm correction)
          -wavelength wl       --  wavelength of the radar wave (m) default=0.236
          -incidence theta     --  incidence angle of the radar wave (degree) default=37
          -range rng           --  range distance from the radar to the center of the interferogram (m) default=866000
          -rms                 --  output RMS of the data misfit grids (mm): rms.grd
          -dem                 --  output DEM error (m): dem.grd

         output:
          disp_##.grd          --  cumulative displacement time series (mm) grids
          vel.grd              --  mean velocity (mm/yr) grids

         example:
          sbas intf.tab scene.tab 88 28 700 1000
        """
        import os
        import subprocess
        import numpy as np
        import math
        import glob
        import datetime

        # cleanup
        for filename in glob.glob(os.path.join(self.basedir, 'disp*.grd')):
            os.remove(filename)
        filename = os.path.join(self.basedir, 'vel.grd')
        if os.path.exists(filename):
            os.remove(filename)

        unwrap = self.open_grids(baseline_pairs[['ref_date', 'rep_date']][:1], 'unwrap')[0]
        dem = self.get_dem(geoloc=True)

        #N=$(wc -l intf.in   | cut -d ' ' -f1)
        #S=$(wc -l scene.tab | cut -d ' ' -f1)

        N = len(baseline_pairs)
        S = len(np.unique(list(baseline_pairs['ref_date']) + list(baseline_pairs['rep_date'])))

        #bounds = self.geoloc().dissolve().envelope.bounds.values[0]
        llmin, ltmin, llmax, ltmax = self.get_master().dissolve().envelope.bounds.values[0].round(3)
        lon0 = (llmin + llmax)/2
        lat0 = (ltmin + ltmax)/2
        elevation0 = float(dem.sel(lat=lat0, lon=lon0, method='nearest'))
        #print ('coords',lon0, lat0, elevation0)
        _,_,_,look_E,look_N,look_U = self.PRM().SAT_look([lon0, lat0, elevation0])
        #print ('satlook', _,_,_,look_E,look_N,look_U)
        incidence = math.atan2(math.sqrt(float(look_E)**2 + float(look_N)**2), float(look_U))*180/np.pi

        ydim, xdim = unwrap.shape

        xmin = int(unwrap.x.min())
        xmax = int(unwrap.x.max())
        near_range, rng_samp_rate, wavelength = self.PRM().get('near_range', 'rng_samp_rate', 'radar_wavelength')
        # calculation below requires bc utility
        rng_pixel_size = 300000000 / rng_samp_rate / 2
        rng = np.round(rng_pixel_size * (xmin+xmax) /2 + near_range)

        intf_tab = self.intftab(baseline_pairs)
        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(intf_tab, 'ascii'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        scene_tab = self.scenetab(baseline_pairs)
        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(scene_tab, 'ascii'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        argv = ['sbas', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}',
                str(N), str(S), str(xdim), str(ydim), '-atm', str(atm), '-smooth', str(smooth),
                '-wavelength', str(wavelength), '-incidence', str(incidence), '-range', str(rng),
                '-rms', '-dem']
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=self.basedir, encoding='ascii')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: sbas', stderr_data)
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: sbas', stdout_data)

        # fix output grid filenames
        for date in np.unique(np.concatenate([baseline_pairs.ref_date,baseline_pairs.rep_date])):
            jdate = datetime.datetime.strptime(date, '%Y-%m-%d').strftime('%Y%j')
            date = date.replace('-','')
            filename1 = os.path.join(self.basedir, f'disp_{jdate}.grd')
            filename2 = os.path.join(self.basedir, f'disp_{date}.grd')
            if os.path.exists(filename1):
                if debug:
                    print ('DEBUG: rename', filename1, filename2)
                os.rename(filename1, filename2)
            #print (jdate, date)

        return
