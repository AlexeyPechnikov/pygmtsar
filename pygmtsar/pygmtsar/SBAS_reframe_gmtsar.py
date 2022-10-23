#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_orbits import SBAS_orbits
from .PRM import PRM

class SBAS_reframe_gmtsar(SBAS_orbits):

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
