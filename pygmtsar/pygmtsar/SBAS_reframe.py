#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_orbits import SBAS_orbits
from .PRM import PRM

class SBAS_reframe(SBAS_orbits):

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

    def get_pins(self, subswath=None):
        """
        Return linear list of two pin coordinates for one or all subswaths. Use this list to easy plot the pins.
        """
    
        if subswath is None:
            # that's ok when pins are not defined here
            # all the pins useful for plotting only
            return self.pins
        else:
            # detect all the subswaths
            subswaths = self.get_subswaths()
            # check pins defined for all the subswaths
            assert len(self.pins)/4. == len(subswaths), f'ERROR: Pins are not defined for all the subswaths. \
                Found {len(self.pins)} pins for {subswaths} subswathss'
            assert subswath in subswaths, f'Subswath {subswath} not found'
            idx = subswaths.tolist().index(subswath)
            pins = self.pins[4*idx:4*(idx+1)]
            assert len(pins) == 4, f'ERROR: wrong number of pins detected. Found {len(pins)} for subswath {subswath}'
            return pins

    def set_pins(self, *args):
        """
        Estimate framed area using two pins on Sentinel-1 GCPs approximation.
        For each defined subswath the pins should be defined like to
            [x1, y1, x2, y2]
            [[x1, y1], [x2, y2]]
            [[x1, y1], None]
            [None, [x2, y2]]
            [None, None]
        The pins automatically ordering for ascending and descending orbits.
        """
        import numpy as np
        from shapely.geometry import Point
        # project pin location to boundary geometry
        def pip2pin(geom, lon, lat):
            pin = geom.boundary.interpolate(geom.boundary.project(Point(lon, lat)))
            return np.round(pin.x, 3), np.round(pin.y, 3)

        # detect all the subswaths
        subswaths = self.get_subswaths()
        if len(args) == 0:
            args = len(subswaths)*[None]
        # check input data to have a single argument for the each subswath
        assert len(args) == len(subswaths), f'Define pair of pins for the each subswath \
            {",".join(map(str,subswaths))} ({2*len(subswaths)} pins and {4*len(subswaths)} coordinates in total)'

        # iterate 1 to 3 subswaths
        allpins = []
        for (subswath, pins) in zip(subswaths, args):
            #print ('subswath, pins', subswath, pins)

            error = False
            warning = False

            if pins is None:
                pins = [None, None]
            if len(pins) == 4:
                pins = np.array(pins).reshape(2,2)
            assert len(pins) == 2, 'Define two pins as two pairs of lat,lon coordinates where pin2 is upper pin1 or use None'
            #print ('pins', pins)
            pin1 = pins[0]
            pin2 = pins[1]

            df = self.get_master(subswath)
            geom = df['geometry'].unary_union
            orbit = df['orbit'][0]

            # check the pins validity
            llmin, ltmin, llmax, ltmax = geom.envelope.bounds
            if not np.all(pin1) is None and not geom.intersects(Point(pin1[0], pin1[1])):
                print (f'ERROR subswath {subswath}: pin1 lays outside of master frame. Move the pin or set it to None and try again.')
                error = True
            if not np.all(pin2) is None and not geom.intersects(Point(pin2[0], pin2[1])):
                print (f'ERROR subswath {subswath}: pin2 lays outside of master frame. Move the pin or set it to None and try again.')
                error = True

            # check pin1
            if np.all(pin1) is None and orbit == 'A':
                # use right bottom corner
                print (f'NOTE subswath {subswath}: pin1 is not defined, master image corner coordinate will be used')
                warning = True
                #pin1 = [llmax, ltmin]
                pin1 = pip2pin(geom, llmax, ltmin)
            elif np.all(pin1) is None and orbit == 'D':
                # use right top corner
                print (f'NOTE subswath {subswath}: pin1 is not defined, master image corner coordinate will be used')
                warning = True
                #pin1 = [llmin, ltmin]
                pin1 = pip2pin(geom, llmin, ltmin)
            # check pin2
            if np.all(pin2) is None and orbit == 'A':
                # use left top corner
                print (f'NOTE subswath {subswath}: pin2 is not defined, master image corner coordinate will be used')
                warning = True
                #pin2 = [llmin, ltmax]
                pin2 = pip2pin(geom, llmin, ltmax)
            elif np.all(pin2) is None and orbit == 'D':
                # use left bottom corner
                print (f'NOTE subswath {subswath}: pin2 is not defined, master image corner coordinate will be used')
                warning = True
                #pin2 = [llmax, ltmax]
                pin2 = pip2pin(geom, llmax, ltmax)

            # check pins order
            if pin1[1] >= pin2[1]:
                print (f'ERROR subswath {subswath}: pin1 is upper than pin2. Fix to set pin1 at bottom and pin2 at top.')
                error = True
            # swap pins for Descending orbit
            if orbit == 'A':
                allpins += [pin1[0], pin1[1], pin2[0], pin2[1]]
            else:
                allpins += [pin2[0], pin2[1], pin1[0], pin1[1]]

        self.pins = allpins
        # pins are defined even is case of errors to have ability to plot them
        assert not error, 'ERROR: Please fix all the issues listed above to continue. Note: you are able to plot the pins.'

    def reframe(self, subswath, date, debug=False):
        """
        Estimate framed area using two pins using Sentinel-1 GCPs approximation.
        """
        import numpy as np
        import shapely
        import os

        pins = self.get_pins(subswath)

        df = self.get_aligned(subswath, date)
        stem = self.multistem_stem(subswath, df['datetime'][0])[1]
        #print ('stem', stem)

        old_filename = os.path.join(self.basedir, f'{stem}')
        #print ('old_filename', old_filename)

        self.make_s1a_tops(subswath, date, debug=debug)

        prm = PRM.from_file(old_filename+'.PRM')
        tmpazi = prm.SAT_llt2rat([pins[0], pins[1], 0], precise=1)[1]
        if debug:
            print ('DEBUG: ','tmpazi', tmpazi)
        prm.shift_atime(tmpazi, inplace=True).update()
        azi1 = prm.SAT_llt2rat([pins[0], pins[1], 0], precise=1)[1] + tmpazi
        azi2 = prm.SAT_llt2rat([pins[2], pins[3], 0], precise=1)[1] + tmpazi
        if debug:
                print ('DEBUG: ','azi1', azi1, 'azi2', azi2)

        # Working on bursts covering $azi1 ($ll1) - $azi2 ($ll2)...
        #print ('assemble_tops', subswath, date, azi1, azi2, debug)
        self.assemble_tops(subswath, date, azi1, azi2, debug=debug)

        # Parse new .xml to define new scene name
        # like to 's1b-iw3-slc-vv-20171117t145922-20171117t145944-008323-00ebab-006'
        filename = os.path.splitext(os.path.split(df['datapath'][0])[-1])[0]
        head1 = filename[:15]
        tail1 = filename[-17:]
        xml_header = self.annotation(old_filename+'.xml')['product']['adsHeader']
        date_new = xml_header['startTime'][:10].replace('-','')
        t1 = xml_header['startTime'][11:19].replace(':','')
        t2 = xml_header['stopTime'][11:19].replace(':','')
        new_name = f'{head1}{date_new}t{t1}-{date_new}t{t2}-{tail1}'
        new_filename = os.path.join(self.basedir, new_name)
        #print ('new_filename', new_filename)

        # rename xml and tiff
        for ext in ['.tiff', '.xml']:
            if debug:
                print('DEBUG: rename', old_filename+ext, new_filename+ext)
            os.rename(old_filename+ext, new_filename+ext)

        # cleanup
        for fname in [old_filename+'.LED', old_filename+'.PRM']:
            if not os.path.exists(fname):
                continue
            if debug:
                print ('DEBUG: remove', fname)
            os.remove(fname)

        # update and return only one record
        df = df.head(1)
        df['datetime'] = self.text2date(f'{date_new}t{t1}', False)
        df['metapath'] = new_filename + '.xml'
        df['datapath'] = new_filename + '.tiff'
        # update approximate location
        gcps = self.geoloc(new_filename + '.xml').geometry
        df['geometry'] = shapely.geometry.MultiPoint(gcps).minimum_rotated_rectangle

        return df

    def reframe_parallel(self, dates=None, n_jobs=-1, **kwargs):
        from tqdm.auto import tqdm
        import joblib
        import pandas as pd

        if len(self.pins) == 0:
            # set auto-pins when the list is empty
            self.set_pins()

        if dates is None:
            dates = self.df.index.unique().values

        subswaths = self.get_subswaths()

        # process all the scenes
        with self.tqdm_joblib(tqdm(desc='Reframing', total=len(dates)*len(subswaths))) as progress_bar:
            records = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.reframe)(subswath, date, **kwargs) \
                                                     for date in dates for subswath in subswaths)

        self.df = pd.concat(records)
