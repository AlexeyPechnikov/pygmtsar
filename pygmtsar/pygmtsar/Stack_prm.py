# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_base import Stack_base
from .PRM import PRM

class Stack_prm(Stack_base):

    def PRM(self, date=None, subswath=None):
        """
        Open a PRM (Parameter) file.

        Parameters
        ----------
        date : str, optional
            The date of the PRM file. If None or equal to self.reference, return the reference PRM file. Default is None.
        multi : bool, optional
            If True, open a multistem PRM file. If False, open a stem PRM file. Default is True.
        singleswath : bool, optional
            If True, open a single-digit subswath PRM file instead of a merged (multi-digit) one. Default is False.

        Returns
        -------
        PRM
            An instance of the PRM class representing the opened PRM file.
        """
        import os

        # check if subswath exists or return a single subswath for None
        if subswath is None:
            subswath = self.get_subswath()

        if date is None:
            date == self.reference

        prefix = self.get_subswath_prefix(subswath, date)
        filename = os.path.join(self.basedir, f'{prefix}.PRM')
        return PRM.from_file(filename)
        
    def PRM_merged(self, date=None, offsets='auto'):
    
        if isinstance(offsets, str) and offsets == 'auto':
            offsets = self.prm_offsets()
            
        maxy, maxx = offsets['extent']
        minh = offsets['bottom']
        return self.PRM(date=date).fix_merged(maxy, maxx, minh)

    def prm_offsets(self, debug=False):
        import xarray as xr
        import numpy as np
        from scipy import constants

        subswaths = self.get_subswaths()
        if not isinstance(subswaths, (str, int)):
            subswaths = ''.join(map(str, subswaths))

        if len(subswaths) == 1:
            prm = self.PRM(subswath=int(subswaths))
            maxx, yvalid, num_patch = prm.get('num_rng_bins', 'num_valid_az', 'num_patches')
            maxy = yvalid * num_patch
            offsets = {'bottom': 0, 'extent': [maxy, maxx]}
            if debug:
                print ('offsets', offsets)
            return offsets

        # calculate the offsets to merge subswaths
        prms = []
        ylims = []
        xlims = []
        for subswath in subswaths:
            #print (subswath)
            prm = self.PRM(subswath=subswath)
            prms.append(prm)
            ylims.append(prm.get('num_valid_az'))
            xlims.append(prm.get('num_rng_bins'))

        assert len(np.unique([prm.get('PRF') for prm in prms])), 'Image PRFs are not consistent'
        assert len(np.unique([prm.get('rng_samp_rate') for prm in prms])), 'Image range sampling rates are not consistent'

        bottoms = [0] + [int(np.round(((prm.get('clock_start') - prms[0].get('clock_start')) * 86400 * prms[0].get('PRF')))) for prm in prms[1:]]
        # head123: 0, 466, -408
        if debug:
            print ('bottoms init', bottoms)
        # minh: -408
        minh = min(bottoms)
        if debug:
            print ('minh', minh)
        #head123: 408, 874, 0
        bottoms = np.asarray(bottoms) - minh
        if debug:
            print ('bottoms', bottoms)

        #ovl12,23: 2690, 2558
        ovls = [prm1.get('num_rng_bins') - \
                int(np.round(((prm2.get('near_range') - prm1.get('near_range')) / (constants.speed_of_light/ prm1.get('rng_samp_rate') / 2)))) \
                for (prm1, prm2) in zip(prms[:-1], prms[1:])]
        if debug:
            print ('ovls', ovls)

        #Writing the grid files..Size(69158x13075)...
        #maxy: 13075
        # for SLC
        maxy = max([prm.get('num_valid_az') + bottom for prm, bottom in zip(prms, bottoms)])
        if debug:
            print ('maxy', maxy)
        maxx = sum([prm.get('num_rng_bins') - ovl - 1 for prm, ovl in zip(prms, [-1] + ovls)])
        if debug:
            print ('maxx', maxx)

        #Stitching location n1 = 1045
        #Stitching location n2 = 935
        ns = [np.ceil(-prm.get('rshift') + prm.get('first_sample') + 150.0).astype(int) for prm in prms[1:]]
        ns = [10 if n < 10 else n for n in ns]
        if debug:
            print ('ns', ns)

        # left and right coordinates for every subswath valid area
        lefts = []
        rights = []

        # 1st
        xlim = prms[0].get('num_rng_bins') - ovls[0] + ns[0]
        lefts.append(0)
        rights.append(xlim)

        # 2nd
        if len(prms) == 2:
            xlim = prms[1].get('num_rng_bins') - 1
        else:
            # for 3 subswaths
            xlim = prms[1].get('num_rng_bins') - ovls[1] + ns[1]
        lefts.append(ns[0])
        rights.append(xlim)

        # 3rd
        if len(prms) == 3:
            xlim = prms[2].get('num_rng_bins') - 2
            lefts.append(ns[1])
            rights.append(xlim)

        # check and merge SLCs
        sumx = sum([right-left for right, left in zip(rights, lefts)])
        if debug:
            print ('assert maxx == sum(...)', maxx, sumx)
        assert maxx == sumx, 'Incorrect output grid range dimension size'

        offsets = {'bottoms': bottoms, 'lefts': lefts, 'rights': rights, 'bottom': minh, 'extent': [maxy, maxx], 'ylims': ylims, 'xlims': xlims}
        if debug:
            print ('offsets', offsets)

        return offsets
