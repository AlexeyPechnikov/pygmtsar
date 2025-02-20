# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_export import Stack_export

class Stack(Stack_export):

    def __init__(self, basedir, drop_if_exists=False):
        """
        Initialize an instance of the Stack class.

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
        Initialize an Stack object with the data directory 'data' and the base directory 'raw':
        stack = Stack('data', basedir='raw')

        Initialize an Stack object with the data directory 'data', DEM filename 'data/DEM_WGS84.nc', and the base directory 'raw':
        stack = Stack('data', 'data/DEM_WGS84.nc', 'raw')
        """
        import os
        os.makedirs(basedir)
        self.basedir = basedir
