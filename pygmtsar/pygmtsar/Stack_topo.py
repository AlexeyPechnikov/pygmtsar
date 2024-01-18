# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_trans_inv import Stack_trans_inv

class Stack_topo(Stack_trans_inv):

    def get_topo(self):
        """
        Get the radar topography grid.

        Returns
        -------
        xarray.DataArray
            The 'topo' variable data is a xarray.DataArray.

        Examples
        --------
        Get DEM:
        topo = stack.get_topo()

        Notes
        -----
        This method returns 'topo' variable from inverse radar transform grid.
        """
        return self.get_trans_inv()['ele'].rename('topo')

    def plot_topo(self, topo='auto', caption='Topography in Radar Coordinates', AOI=None, POI=None, dpi=300, aspect=None):
        import matplotlib.pyplot as plt

        plt.figure(figsize=(12,4), dpi=dpi)
        if isinstance(topo, str) and topo == 'auto':
            self.get_topo().plot.imshow(cmap='gray')
        else:
            topo.plot.imshow(cmap='gray')
        if AOI is not None:
            boundaries = AOI.boundary
            AOI[~boundaries.is_empty].boundary.plot(ax=plt.gca(), color='red')
            AOI[boundaries.is_empty].plot(ax=plt.gca(), color='red')
        if POI is not None:
            POI.plot(ax=plt.gca(), marker='*', markersize=150, color='red')
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        plt.xlabel('Range', fontsize=16)
        plt.ylabel('Azimuth', fontsize=16)
        plt.title(caption, fontsize=18)
        plt.show()

