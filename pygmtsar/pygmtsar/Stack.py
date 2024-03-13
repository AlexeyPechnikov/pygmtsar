# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_export import Stack_export
from .S1 import S1
from .PRM import PRM

class Stack(Stack_export):

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

    @staticmethod
    def plot_AOI(geometry=None, **kwargs):
        import matplotlib.pyplot as plt
        if 'AOI' not in kwargs:
            return
        geometry = kwargs['AOI']
        if geometry is None:
            return
        if 'boundary_color' not in kwargs:
            boundary_color = 'red'
        else:
            boundary_color = kwargs['boundary_color']
        boundaries = geometry.boundary
        geometry[~boundaries.is_empty].boundary.plot(ax=plt.gca(), color=boundary_color)
        geometry[boundaries.is_empty].plot(ax=plt.gca(), color=boundary_color)

    @staticmethod
    def plot_POI(geometry=None, **kwargs):
        import matplotlib.pyplot as plt
        if 'POI' not in kwargs:
            return
        geometry = kwargs['POI']
        if geometry is None:
            return
        if 'marker' not in kwargs:
            marker = '*'
        else:
            marker = kwargs['marker']
        if 'marker_size' not in kwargs:
            marker_size = 100
        else:
            marker_size = kwargs['marker_size']
        if 'marker_color' not in kwargs:
            marker_color = 'red'
        else:
            marker_color = kwargs['marker_color']
        geometry.plot(ax=plt.gca(), marker=marker, markersize=marker_size, color=marker_color)

    def plot_scenes(self, dem='auto', caption='Estimated Scene Locations', cmap='turbo', aspect=None, **kwargs):
        import matplotlib.pyplot as plt
        import matplotlib

        plt.figure()
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
        # too small an alpha becomes invisible
        alpha = 0.5/len(gdf)
        alpha = alpha if alpha>=0.002 else 0.002
        gdf.reset_index().plot(color=[colors[k] for k in gdf.index], alpha=alpha, edgecolor='black', ax=plt.gca())
        self.plot_AOI(**kwargs)
        self.plot_POI(**kwargs)
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        plt.title(caption)
