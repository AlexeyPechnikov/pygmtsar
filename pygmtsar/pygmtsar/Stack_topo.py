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

    def plot_topo(self, topo='auto', caption='Topography on WGS84 ellipsoid, [m]',
                  quantile=None, vmin=None, vmax=None, cmap='gray', aspect=None, **kwargs):
        import numpy as np
        import matplotlib.pyplot as plt

        if isinstance(topo, str) and topo == 'auto':
            topo = self.get_topo()
            
        if quantile is not None:
            assert vmin is None and vmax is None, "ERROR: arguments 'quantile' and 'vmin', 'vmax' cannot be used together"

        if quantile is not None:
            vmin, vmax = np.nanquantile(topo, quantile)

        plt.figure()
        topo.plot.imshow(cmap=cmap, vmin=vmin, vmax=vmax)
        self.plot_AOI(**kwargs)
        self.plot_POI(**kwargs)
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        plt.xlabel('Range')
        plt.ylabel('Azimuth')
        plt.title(caption)

