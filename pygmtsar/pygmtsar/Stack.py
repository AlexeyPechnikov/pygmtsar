# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_ps import Stack_ps
from .S1 import S1
from .PRM import PRM

class Stack(Stack_ps):

    df = None
    basedir = None
    reference = None
    dem_filename = None
    landmask_filename = None
    
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
        import shutil

        # (re)create basedir only when force=True
        if os.path.exists(basedir):
            if drop_if_exists:
                shutil.rmtree(basedir)
            else:
                raise ValueError('ERROR: The base directory already exists. Use drop_if_exists=True to delete it and start new processing.')
        os.makedirs(basedir)
        self.basedir = basedir

    def set_scenes(self, scenes):
        assert len(scenes), 'ERROR: the scenes list is empty.'
        assert len(scenes[scenes.orbitpath.isna()])==0, 'ERROR: orbits missed, check "orbitpath" column.'
        self.df = scenes
        if self.reference is None:
            print (f'NOTE: auto set reference scene {scenes.index[0]}. You can change it like Stack.set_reference("2022-01-20")')
            self.reference = self.df.index[0]
        return self

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

    def plot_scenes(self, AOI=None, POI=None, dem='auto', caption='Estimated Scene Locations', cmap='turbo', dpi=150, aspect=None):
        import matplotlib.pyplot as plt
        import matplotlib

        plt.figure(figsize=(12, 4), dpi=dpi)
        if isinstance(dem, str) and dem == 'auto':
            if self.dem_filename is not None:
                dem = self.get_dem()
                # TODO: check shape and decimate large grids
                dem.plot.imshow(cmap='gray', add_colorbar=False)
        elif dem is not None:
            dem.plot.imshow(cmap='gray', add_colorbar=False)
        gdf = self.to_dataframe()
        cmap = matplotlib.colormaps[cmap]
        colors = dict([(v, cmap(k)) for k, v in enumerate(gdf.index.unique())])
        gdf.reset_index().plot(color=[colors[k] for k in gdf.index], alpha=0.5/len(gdf), edgecolor='black', ax=plt.gca())
        if AOI is not None:
            boundaries = AOI.boundary
            AOI[~boundaries.is_empty].boundary.plot(ax=plt.gca(), color='red')
            AOI[boundaries.is_empty].plot(ax=plt.gca(), color='red')
        if POI is not None:
            POI.plot(ax=plt.gca(), marker='*', markersize=150, color='red')
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        plt.title(caption, fontsize=18)
        plt.show()
