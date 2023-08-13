# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_ps import SBAS_ps
from .S1 import S1
from .PRM import PRM

class SBAS(SBAS_ps):

    def __init__(self, basedir, scenes, reference=None, dem_filename=None, landmask_filename=None, drop_if_exists=False):
        """
        Initialize an instance of the SBAS class.

        Parameters
        ----------
        basedir : str
            The base directory for processing.
        scenes : GeoPandas Dataframe
            Sentinel-1 scenes with bursts geometries and orbits in structured format. 
        dem_filename : str, optional
            The filename of the DEM (Digital Elevation Model) WGS84 NetCDF file. Default is None.
        landmask_filename : str, optional
            The filename of the landmask WGS84 NetCDF file. Default is None.

        Examples
        --------
        Initialize an SBAS object with the data directory 'data' and the base directory 'raw':
        sbas = SBAS('data', basedir='raw')

        Initialize an SBAS object with the data directory 'data', DEM filename 'data/DEM_WGS84.nc', and the base directory 'raw':
        sbas = SBAS('data', 'data/DEM_WGS84.nc', 'raw')
        """
        import os
        import shutil

        # (re)create basedir only when force=True
        if os.path.exists(basedir):
            if drop_if_exists:
                shutil.rmtree(basedir)
            else:
                raise ValueError('ERROR: The base directory already exists. Use drop_if_exists=True to delete it and start new processing.')
        os.makedirs(basedir)
        self.basedir = basedir

        self.df = scenes
        if reference is None:
            print (f'NOTE: reference scene is not defined, use {scenes.index[0]}. You can change it like SBAS.set_reference("2022-01-20")')
            self.reference = self.df.index[0]
        self.set_dem(dem_filename)
        self.set_landmask(landmask_filename)

#    def make_gaussian_filter(self, range_dec, azi_dec, wavelength, debug=False):
#        """
#        Wrapper for PRM.make_gaussian_filter() and sonamed command line tool. Added for development purposes only.
#        """
#        import numpy as np
#
#        gauss_dec, gauss_string = self.PRM().make_gaussian_filter(range_dec, azi_dec, wavelength, debug=debug)
#        coeffs = [item for line in gauss_string.split('\n') for item in line.split('\t') if item != '']
#        # x,y dims order
#        shape = np.array(coeffs[0].split(' ')).astype(int)
#        # y,x dims order
#        matrix = np.array(coeffs[1:]).astype(float).reshape((shape[1],shape[0]))
#        return (gauss_dec, matrix)
