# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_multilooking import Stack_multilooking

class Stack_landmask(Stack_multilooking):

    def set_landmask(self, landmask_filename):
        """
        Set the landmask file path.

        Parameters
        ----------
        landmask_filename : str or None
            The NetCDF file path of the landmask. If None, the landmask file path will be set to None.

        Examples
        --------
        stack = stack.set_landmask('data/landmask.nc')

        Returns
        -------
        self : Stack
            The Stack object with the updated landmask file path.
        """
        import os
        if landmask_filename is not None:
            assert os.path.exists(landmask_filename), f'Landmask file not found: {landmask_filename}'
            assert os.path.isfile(landmask_filename) and os.access(landmask_filename, os.R_OK), f'Landmask file is not readable: {landmask_filename}'
            self.landmask_filename = os.path.relpath(landmask_filename,'.')
        else:
            self.landmask_filename = None
        return self

    def get_landmask(self):
        """
        Get the landmask grid in geographic coordinates.

        Parameters
        ----------
        None

        Returns
        -------
        landmask : xarray.DataArray
            The landmask grid.

        Examples
        --------
        Get land mask in geographic coordinates:
        landmask = stack.get_landmask()

        Notes
        -----
        This method opens the landmask NetCDF file and extracts the landmask grid.
        """
        import xarray as xr
        import os

        if self.landmask_filename is None:
            raise Exception('Set landmask first')

        # open DEM file and find the elevation variable
        # because sometimes grid includes 'crs' or other variables
        landmask = xr.open_dataset(self.landmask_filename, engine=self.netcdf_engine, chunks=self.chunksize)
        assert 'lat' in landmask.coords and 'lon' in landmask.coords
        # define latlon array
        z_array_name = [data_var for data_var in landmask.data_vars if len(landmask.data_vars[data_var].coords)==2]
        assert len(z_array_name) == 1
        # extract the array and fill missed values by zero (mostly ocean area)
        landmask = landmask[z_array_name[0]].fillna(0)
        # round the coordinates up to 1 mm to have the same grid for dem and landmask
        landmask['lat'] = landmask.lat.round(8)
        landmask['lon'] = landmask.lon.round(8)
        return landmask

    def download_landmask(self, product='1s', debug=False):
        """
        Download the landmask and save as NetCDF file.

        Parameters
        ----------
        product : str, optional
                Available options are '1s' (1 arcsec ~= 30m, default) and '3s' (3 arcsec ~= 90m).
        debug : bool, optional
            Whether to enable debug mode. Default is False.

        Examples
        --------
        stack.download_landmask()

        Notes
        -----
        This method downloads the landmask using GMT's local data or server. The landmask is built based on the interferogram DEM area.
        """
        import pygmt
        import os
        from tqdm.auto import tqdm

        if self.landmask_filename is not None:
            print ('NOTE: landmask exists, ignore the command. Use Stack.set_landmask(None) to allow new landmask downloading')
            return

        # generate the same as DEM grid
        landmask_filename = os.path.join(self.basedir, 'landmask.nc')

        # geographical grid for interferogram area only
        dem = self.get_dem()
        llmin = dem.lon.min().item()
        llmax = dem.lon.max().item()
        ltmin = dem.lat.min().item()
        ltmax = dem.lat.max().item()
        region = f'{llmin}/{llmax}/{ltmin}/{ltmax}'
        if debug:
            print('region', region)

        with tqdm(desc='Landmask Downloading', total=1) as pbar:
            landmask = pygmt.grdlandmask(resolution='f', region=region, spacing=product, maskvalues='NaN/1')
            if os.path.exists(landmask_filename):
                os.remove(landmask_filename)
            landmask.to_netcdf(landmask_filename)
            pbar.update(1)

        self.landmask_filename = landmask_filename

    def plot_landmask(self, landmask='auto', caption='Land Mask', cmap='binary', aspect=None, **kwargs):
        import matplotlib.pyplot as plt

        if isinstance(landmask, str) and landmask == 'auto':
            landmask = self.get_landmask()

        plt.figure()
        landmask.plot.imshow(vmin=0, cmap=cmap)
        self.plot_AOI(**kwargs)
        self.plot_POI(**kwargs)
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        plt.title(caption)
