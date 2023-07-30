# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_trans_inv import SBAS_trans_inv
from .PRM import PRM
from .tqdm_dask import tqdm_dask

class SBAS_topo_ra(SBAS_trans_inv):

    def topo_ra(self, subswath=None, chunksize=None, interactive=False):
        """
        Compute the topography in radar coordinates (topo_ra).

        Parameters
        ----------
        subswath : int or None, optional
            The subswath number to compute the topographic radar coordinates for. If None, the computation
            will be performed for all subswaths. Default is None.
        interactive : bool, optional
            If True, the computation will be performed interactively and the result will be returned as a delayed object.
            Default is False.
        """
        import dask
        import xarray as xr
        import numpy as np
        import os

        # GMTSAR phasediff tool requires "The dimension SLC must be multiplication factor of the topo_ra"
        topo_ra = self.get_trans_dat_inv(subswath).ele[1:,1:].rename('topo_ra')
    
        if interactive:
            # do not flip vertically because it's returned as is without SBAS.get_topo_ra() function
            return topo_ra

        # save to NetCDF file
        filename = self.get_filenames(subswath, None, 'topo_ra')
        if os.path.exists(filename):
            os.remove(filename)
        # flip vertically for GMTSAR compatibility reasons
        topo_ra = xr.DataArray(dask.array.flipud(topo_ra), coords=topo_ra.coords, name=topo_ra.name)
        # rename to save lazy NetCDF preventing broken coordinates (y,y) 
        topo_ra = topo_ra.rename({'y': 'a', 'x': 'r'})
        handler = topo_ra.to_netcdf(filename,
                                    encoding={'topo_ra': self.compression(topo_ra.shape, chunksize=chunksize)},
                                    engine=self.engine,
                                    compute=False)
        return handler

    def topo_ra_parallel(self, coarsen=None, interactive=False, **kwargs):
        """
        Build topography in radar coordinates from WGS84 DEM using parallel computation.

        Parameters
        ----------
        interactive : bool, optional
            If True, the computation will be performed interactively and the results will be returned as delayed objects.
            If False, the progress will be displayed using tqdm_dask. Default is False.

        Returns
        -------
        handler or list of handlers
            The handler(s) of the delayed computation if 'interactive' is True. Otherwise, None.

        Examples
        --------
        sbas.topo_ra_parallel()

        Notes
        -----
        This method performs the parallel computation of topography in the radar coordinates (topo_ra) for all subswaths
        using Dask. It calls the 'topo_ra' method for each subswath in parallel. If 'interactive' is True, the delayed
        computation handlers will be returned. Otherwise, the progress will be displayed using tqdm_dask.
        """
        import dask

        # default argument value None is required to call it from geocode_parallel function
        if coarsen is None:
            coarsen = 4
            print (f'Note: use default transform grid spacing {coarsen}')

        # generate the coordinates transform
        self.trans_dat_parallel(coarsen=coarsen, **kwargs)
        # generate the inversion transform
        self.trans_dat_inv_parallel(coarsen=coarsen, **kwargs)

        # process all the subswaths
        subswaths = self.get_subswaths()
        delayeds = []
        for subswath in subswaths:
            delayed = self.topo_ra(subswath=subswath, interactive=interactive, **kwargs)
            if not interactive:
                tqdm_dask(dask.persist(delayed), desc=f'Radar Topography Computing sw{subswath}')
            else:
                delayeds.append(delayed)

        if interactive:
            return delayeds[0] if len(delayeds)==1 else delayeds

        # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
        import gc; gc.collect()

    def get_topo_ra(self, chunksize=None):
        """
        Get the radar topography grid.

        Returns
        -------
        xarray.DataArray or list of xarray.DataArray
            The 'topo_ra' grid data as a single xarray.DataArray if only one grid is found,
            or a list of xarray.DataArray if multiple grids are found.

        Examples
        --------
        Get DEM for all the processed subswaths:
        topo_ra = sbas.get_topo_ra()

        Notes
        -----
        This method opens the 'topo_ra' grids using the `open_grids` method and applies the `func` function
        to each grid to flip it vertically for compatibility reasons with GMTSAR. The resulting grids are returned.
        """
        import xarray as xr
        import dask.array

        if chunksize is None:
            chunksize = self.chunksize

        def func(topo):
            # flip vertically for GMTSAR compatibility reasons
            topo = xr.DataArray(dask.array.flipud(topo), coords=topo.coords, attrs=topo.attrs, name=topo.name)
            # fix renamed dimensions required to save lazy NetCDF properly
            if 'a' in topo.dims and 'r' in topo.dims:
                return topo.rename({'a': 'y', 'r': 'x'})
            return topo

        topos = self.open_grids(None, 'topo_ra', chunksize=chunksize, func=func)

        return topos[0] if len(topos)==1 else topos
