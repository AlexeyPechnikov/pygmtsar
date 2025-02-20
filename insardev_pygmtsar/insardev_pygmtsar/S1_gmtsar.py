# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .S1_prm import S1_prm
from .PRM import PRM

class S1_gmtsar(S1_prm):

    def _ext_orb_s1a(self, burst, date=None, debug=False):
        """
        Extracts orbital data for the Sentinel-1A satellite by running GMTSAR binary `ext_orb_s1a`.

        Parameters
        ----------
        stem : str
            Stem name used for file naming.
        date : str, optional
            Date for which to extract the orbital data. If not provided or if date is the reference, 
            it will extract the orbital data for the reference. Defaults to None.
        debug : bool, optional
            If True, prints debug information. Defaults to False.

        Examples
        --------
        _ext_orb_s1a(1, 'stem_name', '2023-05-24', True)
        """
        import os
        import subprocess

        if date is None or date == self.reference:
            date == self.reference
            df = self.get_reference(burst)
        else:
            df = self.get_repeat(burst, date)
        #print ('df', df)

        prefix = self.get_prefix(burst, date)
        #print ('prefix', prefix)
        if os.path.dirname(prefix) == '':
            basedir = self.basedir
        else:
            basedir = os.path.join(self.basedir, os.path.dirname(prefix))
            prefix = os.path.basename(prefix)
        #print ('basedir', basedir)

        path = df['path'].iloc[0]
        orbit = df['orbit'].iloc[0]
        orbitfile = os.path.join(path, orbit)
        orbitfile = os.path.relpath(orbitfile, basedir)

        argv = ['ext_orb_s1a', f'{prefix}.PRM', orbitfile, prefix]
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8', cwd=basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: ext_orb_s1a', stderr_data)
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: ext_orb_s1a', stdout_data)

        return

    # produce LED and PRM in basedir
    # when date=None work on reference scene
    def _make_s1a_tops(self, burst, date=None, mode=0, rshift_fromfile=None, ashift_fromfile=None, debug=False):
        """
        Produces LED and PRM in the base directory by executing GMTSAR binary `make_s1a_tops`.

        Parameters
        ----------
        date : str, optional
            Date for which to create the Sentinel-1A TOPS products. If not provided, 
            it processes the reference image. Defaults to None.
        mode : int, optional
            Mode for `make_s1a_tops` script: 
            0 - no SLC; 
            1 - center SLC; 
            2 - high SLCH and low SLCL; 
            3 - output ramp phase.
            Defaults to 0.
        rshift_fromfile : str, optional
            Path to the file with range shift data. Defaults to None.
        ashift_fromfile : str, optional
            Path to the file with azimuth shift data. Defaults to None.
        debug : bool, optional
            If True, prints debug information. Defaults to False.

        Notes
        -----
        The function executes an external binary `make_s1a_tops`.
        Also, this function calls the `ext_orb_s1a` method internally.

        Examples
        --------
        _make_s1a_tops(1, '2023-05-24', 1, '/path/to/rshift.grd', '/path/to/ashift.grd', True)
        """
        import os
        import subprocess

        #or date == self.reference
        if date is None:
            date = self.reference
            df = self.get_reference(burst)
            # for reference image mode should be 1
            mode = 1
        else:
            df = self.get_repeat(burst, date)
        
        prefix = self.get_prefix(burst, date)
        if os.path.dirname(prefix) == '':
            basedir = self.basedir
        else:
            basedir = os.path.join(self.basedir, os.path.dirname(prefix))
            prefix = os.path.basename(prefix)
            if rshift_fromfile is not None:
                rshift_fromfile = os.path.basename(rshift_fromfile)
            if ashift_fromfile is not None:
                ashift_fromfile = os.path.basename(ashift_fromfile)

        name = df['burst'].iloc[0]
        path = df['path'].iloc[0]
        xmlfile = os.path.join(path, burst, 'annotation', name + '.xml')
        xmlfile = os.path.relpath(xmlfile, basedir)
        tiffile = os.path.join(path, burst, 'measurement', name + '.tiff')
        tiffile = os.path.relpath(tiffile, basedir)

        argv = ['make_s1a_tops', xmlfile, tiffile, prefix, str(mode)]
        if rshift_fromfile is not None:
            argv.append(rshift_fromfile)
        if ashift_fromfile is not None:
            argv.append(ashift_fromfile)
        if debug:
            print ('DEBUG: argv', argv)
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8', cwd=basedir)
        stdout_data, stderr_data = p.communicate()
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: make_s1a_tops', stderr_data)
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: make_s1a_tops', stdout_data)

        self._ext_orb_s1a(burst, date, debug=debug)

        return

