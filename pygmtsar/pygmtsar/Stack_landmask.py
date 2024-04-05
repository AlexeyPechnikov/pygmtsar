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
        print ('NOTE: Function is removed. Download land mask using Tiles().download_landmask()')
        print ('and load with Stack.load_landmask() function.')

    def load_landmask(self, data, geometry='auto'):
        """
        Load and preprocess land mask data from specified datafile or variable.

        Parameters
        ----------
        data : xarray dataarray or str
            Land mask filename or variable.

        Returns
        -------
        None

        Examples
        --------
        Load and crop from local NetCDF file:
        stack.load_landmask('landmask.nc')

        Load and crop from local GeoTIF file:
        stack.load_landmask('landmask.tif')

        Load from Xarray DataArray or Dataset:
        stack.load_landmask(None).load_landmask(landmask)
        stack.load_landmask(None).load_landmask(landmask.to_dataset())
        
        """
        import xarray as xr
        import numpy as np
        import rioxarray as rio
        import geopandas as gpd
        import os

        # generate the same as DEM grid
        landmask_filename = os.path.join(self.basedir, 'landmask.nc')
        
        if self.landmask_filename is not None:
            print ('NOTE: landmask exists, ignore the command. Use Stack.set_landmask(None) to allow new landmask downloading')
            return

        if isinstance(data, (xr.Dataset)):
            landmask = data[list(data.data_vars)[0]]
        elif isinstance(data, (xr.DataArray)):
            landmask = data
        elif isinstance(data, str) and os.path.splitext(data)[-1] in ['.tiff', '.tif', '.TIF']:
            landmask = rio.open_rasterio(data, chunks=self.chunksize).squeeze(drop=True)\
                .rename({'y': 'lat', 'x': 'lon'})\
                .drop('spatial_ref')
            if landmask.lat.diff('lat')[0].item() < 0:
                landmask = landmask.reindex(lat=landmask.lat[::-1])
        elif isinstance(data, str) and os.path.splitext(data)[-1] in ['.nc', '.netcdf', '.grd']:
            landmask = xr.open_dataarray(data, engine=self.netcdf_engine, chunks=self.chunksize)
        elif isinstance(data, str):
            print ('ERROR: filename extension is not recognized. Should be one from .tiff, .tif, .TIF, .nc, .netcdf, .grd')
        else:
            print ('ERROR: argument is not an Xarray object and it is not a file name')

        # unify to DEM
        dem = self.get_dem()
        landmask = landmask.transpose('lat','lon').reindex_like(dem, method='nearest')

        if os.path.exists(landmask_filename):
            os.remove(landmask_filename)
        encoding = {'landmask': self._compression(landmask.shape)}
        landmask.rename('landmask').load()\
            .to_netcdf(landmask_filename, encoding=encoding, engine=self.netcdf_engine)

        self.landmask_filename = landmask_filename

    def plot_landmask(self, landmask='auto', caption='Land Mask', cmap='binary_r', aspect=None, **kwargs):
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
