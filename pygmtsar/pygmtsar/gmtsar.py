#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar

class gmtsar:

    def assemble_tops(self, subswath, date, azi_1, azi_2, debug=False):
        """
        Usage: assemble_tops azi_1 azi_2 name_stem1 name_stem2 ... output_stem

        Example: assemble_tops 1685 9732 s1a-iw1-slc-vv-20150706t135900-20150706t135925-006691-008f28-001
            s1a-iw1-slc-vv-20150706t135925-20150706t135950-006691-008f28-001
            s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001

        Output:s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001.xml
            s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001.tiff

        Note: output files are bursts that covers area between azi_1 and azi_2, set them to 0s to output all bursts
        """
        import numpy as np
        import os
        import subprocess

        df = self.get_aligned(subswath, date)
        #print ('scenes', len(df))

        # assemble_tops requires the same path to xml and tiff files
        datadirs = [os.path.split(path)[:-1] for path in df['datapath']]
        metadirs = [os.path.split(path)[:-1] for path in df['metapath']]
        if not datadirs == metadirs:
            # in case when the files placed in different directories we need to create symlinks for them
            datapaths = []
            for datapath, metapath in zip(df['datapath'], df['metapath']):
                for filepath in [datapath, metapath]:
                    filename = os.path.split(filepath)[-1]
                    relname = os.path.join(self.basedir, filename)
                    if os.path.exists(relname) or os.path.islink(relname):
                        os.remove(relname)
                    os.symlink(os.path.relpath(filepath, self.basedir), relname)
                datapaths.append(os.path.splitext(filename)[0])
        else:
            datapaths = [os.path.relpath(path, self.basedir)[:-5] for path in df['datapath']]
        #print ('datapaths', datapaths)
        stem = self.multistem_stem(subswath, df['datetime'][0])[1]

        # round values and convert to strings
        azi_1 = np.round(azi_1).astype(int).astype(str)
        azi_2 = np.round(azi_2).astype(int).astype(str)

        argv = ['assemble_tops', azi_1, azi_2] + datapaths + [stem]
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: assemble_tops', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: assemble_tops', stdout_data.decode('ascii'))

        return

    def ext_orb_s1a(self, subswath, stem, date=None, debug=False):
        import os
        import subprocess

        if date is None or date == self.master:
            df = self.get_master(subswath)
        else:
            df = self.get_aligned(subswath, date)

        orbit = os.path.relpath(df['orbitpath'][0], self.basedir)

        argv = ['ext_orb_s1a', f'{stem}.PRM', orbit, stem]
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: ext_orb_s1a', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: ext_orb_s1a', stdout_data.decode('ascii'))

        return
    
    # produce LED and PRM in basedir
    # when date=None work on master image
    def make_s1a_tops(self, subswath, date=None, mode=0, rshift_fromfile=None, ashift_fromfile=None, debug=False):
        """
        Usage: make_slc_s1a_tops xml_file tiff_file output mode dr.grd da.grd
         xml_file    - name of xml file 
         tiff_file   - name of tiff file 
         output      - stem name of output files .PRM, .LED, .SLC 
         mode        - (0) no SLC; (1) center SLC; (2) high SLCH and lowSLCL; (3) output ramp phase
        """
        import os
        import subprocess

        #or date == self.master
        if date is None:
            df = self.get_master(subswath)
            # for master image mode should be 1
            mode = 1
        else:
            df = self.get_aligned(subswath, date)

        # TODO: use subswath
        xmlfile = os.path.relpath(df['metapath'][0], self.basedir)
        datafile = os.path.relpath(df['datapath'][0], self.basedir)
        stem = self.multistem_stem(subswath, df['datetime'][0])[1]

        argv = ['make_s1a_tops', xmlfile, datafile, stem, str(mode)]
        if rshift_fromfile is not None:
            argv.append(rshift_fromfile)
        if ashift_fromfile is not None:
            argv.append(ashift_fromfile)
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: make_s1a_tops', stderr_data.decode('ascii'))
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: make_s1a_tops', stdout_data.decode('ascii'))

        self.ext_orb_s1a(subswath, stem, date, debug=debug)

        return

    # stem_tofile + '.PRM' generating
    def merge_swath(self, conf, grid_tofile, stem_tofile, debug=False):
        import subprocess

        argv = ['merge_swath', '/dev/stdin', grid_tofile, stem_tofile]
        if debug:
            print ('DEBUG: argv', argv)
            print ('DEBUG: conf:', f'\n{conf}')
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='ascii')
        stdout_data, stderr_data = p.communicate(input=conf)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: merge_swath', stderr_data)
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: merge_swath', stdout_data)

        return

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
        prm = self.PRM()

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
        _,_,_,look_E,look_N,look_U = prm.SAT_look([lon0, lat0, elevation0])
        #print ('satlook', _,_,_,look_E,look_N,look_U)
        incidence = math.atan2(math.sqrt(float(look_E)**2 + float(look_N)**2), float(look_U))*180/np.pi

        ydim, xdim = unwrap.shape

        xmin = int(unwrap.x.min())
        xmax = int(unwrap.x.max())
        near_range, rng_samp_rate, wavelength = prm.get('near_range', 'rng_samp_rate', 'radar_wavelength')
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
