# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_phasediff import Stack_phasediff

class Stack_landmask(Stack_phasediff):

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
            self.landmask_filename = os.path.relpath(landmask_filename,'.')
        else:
            self.landmask_filename = None
        return self

    def get_landmask(self, inverse_geocode=False):
        """
        Get the landmask grid in geographic or radar coordinates.

        Parameters
        ----------
        inverse_geocode : bool, optional
            Whether to perform inverse geocoding on the landmask grid. If True, the grid will be transformed from geographic
            coordinates to radar coordinates. Default is False.

        Returns
        -------
        landmask : xarray.DataArray
            The landmask grid.

        Examples
        --------
        Get land mask in geographic coordinates:
        landmask_ll = stack.get_landmask()

        Get land mask in radar coordinates:
        landmask_ra = stack.get_landmask(inverse_geocode=True)

        Notes
        -----
        This method opens the landmask NetCDF file and extracts the landmask grid. The grid can be cropped to the valid area
        only or transformed from geographic coordinates to radar coordinates using inverse geocoding.
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

        # double conversion to unify landmask to interferogram grid in radar coordinates
        # inverse geocode to radar coordinates grid
        landmask_ra = self.intf_ll2ra(landmask)
        if inverse_geocode:
            return landmask_ra

        # geocode to the exact geographical coordinates grid as interferogram
        landmask_ll = self.intf_ra2ll(landmask_ra)
        return landmask_ll

    def download_landmask(self, backend=None, debug=False):
        """
        Download the landmask and save as NetCDF file.

        Parameters
        ----------
        backend : deprecated
            This parameter is deprecated and should be omitted.
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
        # correspond to SRTM3 DEM
        arcsec_degree = 0.000833333333333/3

        if self.landmask_filename is not None:
            print ('NOTE: landmask exists, ignore the command. Use Stack.set_landmask(None) to allow new landmask downloading')
            return

        if backend is not None:
            print ('Note: backend argument is deprecated, just omit it')

        # generate the same as DEM grid
        landmask_filename = os.path.join(self.basedir, 'landmask.nc')

        # geographical grid for interferogram area only
        dem = self.get_intf_ra2ll()[['lat','lon']]
        scale = dem.lon.diff('lon')[0].item()
        llmin = dem.lon.min().item()
        llmax = dem.lon.max().item()
        ltmin = dem.lat.min().item()
        ltmax = dem.lat.max().item()
        region = f'{llmin}/{llmax}/{ltmin}/{ltmax}'
        #print('region', region)

        # define approximate resolution in arc seconds
        spacing = (dem.lat.diff('lat')[0]/arcsec_degree).round(3).item()
        # format as string
        spacing = f'{spacing}s'

        with tqdm(desc='Landmask Downloading', total=1) as pbar:
            landmask = pygmt.grdlandmask(resolution='f', region=region, spacing=spacing, maskvalues='NaN/1')
            if os.path.exists(landmask_filename):
                os.remove(landmask_filename)
            landmask.to_netcdf(landmask_filename)
            pbar.update(1)

        self.landmask_filename = landmask_filename
