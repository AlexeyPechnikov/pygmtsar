# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .S1_slc import S1_slc
from .PRM import PRM

class S1_prm(S1_slc):

    def PRM(self, burst, date=None):
        """
        Open a PRM (Parameter) file.

        Parameters
        ----------
        date : str, optional
            The date of the PRM file. If None or equal to self.reference, return the reference PRM file. Default is None.
        multi : bool, optional
            If True, open a multistem PRM file. If False, open a stem PRM file. Default is True.
        
        Returns
        -------
        PRM
            An instance of the PRM class representing the opened PRM file.
        """
        import os

        # TODO
        assert len(burst)!=10, 'ERROR: mixed burst and date arguments (burst={burst} date={date})'

        if date is None:
            date == self.reference

        prefix = self.get_prefix(burst, date)
        filename = os.path.join(self.basedir, f'{prefix}.PRM')
        #print ('PRM filename', filename)
        return PRM.from_file(filename)

    def prm_offsets(self, burst, debug=False):
        import xarray as xr
        import numpy as np
        from scipy import constants

        prm = self.PRM(burst)
        maxx, yvalid, num_patch = prm.get('num_rng_bins', 'num_valid_az', 'num_patches')
        maxy = yvalid * num_patch
        offsets = {'extent': [maxy, maxx]}
        if debug:
            print ('offsets', offsets)
        return offsets

