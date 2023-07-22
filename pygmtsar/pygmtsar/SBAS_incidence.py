# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_geocode import SBAS_geocode
from .tqdm_dask import tqdm_dask

class SBAS_incidence(SBAS_geocode):

    def get_sat_look(self):
        """
        Return satellite look vectors in geographic coordinates as Xarray Dataset.

        Returns
        -------
        xarray.Dataset
            The satellite look vectors in geographic coordinates.

        Examples
        --------
        Get satellite look vectors:
        sat_look_ll = sbas.get_sat_look()

        Notes
        -----
        This function returns the satellite look vectors in geographic coordinates as Xarray Dataset. The satellite look vectors
        should be computed and saved prior to calling this function using the `sat_look_parallel` method.
        """
        import xarray as xr
        import os

        sat_look_file = self.get_filenames(None, None, 'sat_look')
        assert os.path.exists(sat_look_file), 'ERROR: satellite looks grid missed. Build it first using SBAS.sat_look_parallel()'
        sat_look = xr.open_dataset(sat_look_file, engine=self.engine, chunks=self.chunksize).rename({'yy': 'lat', 'xx': 'lon'})

        return sat_look

    #gmt grdmath unwrap_mask.grd $wavel MUL -79.58 MUL = los.grd
    def los_displacement_mm(self, unwraps):
        """
        Compute line-of-sight (LOS) displacement in millimeters.

        Parameters
        ----------
        unwraps : xarray.DataArray or xarray.Dataset
            Unwrapped phase grid(s) in radar or geographic coordinates.

        Returns
        -------
        xarray.DataArray
            Line-of-sight (LOS) displacement grid(s) in millimeters.

        Examples
        --------
        Calculate LOS displacement for unwrapped phase grids in radar coordinates:
        unwraps_ra = sbas.open_grids(pairs, 'unwrap')
        los_disp_ra = sbas.los_displacement_mm(unwraps_ra)
        # or the same code in one line
        los_disp_ra = sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.

        Calculate LOS displacement for detrended unwrapped phase grids in geographic coordinates:
        detrend_ll = sbas.open_grids(pairs, 'detrend', geocode=True)
        los_disp_ll = sbas.los_displacement_mm(detrend_ll)
        # or the same code in one line
        los_disp_ll = sbas.open_grids(pairs, 'detrend', geocode=True, func=sbas.los_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.
        """
        # constant is negative to make LOS = -1 * range change
        # constant is (1000 mm) / (4 * pi)
        scale = -79.58 * self.PRM().get('radar_wavelength')
        los_disp = scale*unwraps
        return los_disp.rename('los')

    def incidence_angle(self):
        """
        Compute the incidence angle grid in geographic coordinates.

        Returns
        -------
        xarray.DataArray
            The incidence angle grid in geographic coordinates.

        Examples
        --------
        Compute the incidence angle grid:
        inc_angle_ll = sbas.incidence_angle()

        Notes
        -----
        This function computes the incidence angle grid in geographic coordinates based on the satellite look vectors.
        The satellite look vectors should be computed and saved prior to calling this function using the `sat_look_parallel` method.
        The incidence angle is calculated using the formula:
        incidence_angle = arctan2(sqrt(look_E**2 + look_N**2), look_U)
        """
        import xarray as xr
        import numpy as np

        sat_look = self.get_sat_look()
        incidence_ll = np.arctan2(np.sqrt(sat_look.look_E**2 + sat_look.look_N**2), sat_look.look_U).rename('incidence_angle')
        return incidence_ll

    def vertical_displacement_mm(self, unwraps):
        """
        Compute vertical displacement in millimeters in geographic coordinates.

        Parameters
        ----------
        unwraps : xarray.DataArray or xarray.Dataset
            Unwrapped phase grid(s) in geographic coordinates.

        Returns
        -------
        xarray.DataArray
            Vertical displacement grid(s) in millimeters.

        Examples
        --------
        Calculate vertical displacement for unwrapped phase grids in geographic coordinates:
        unwraps_ll = sbas.open_grids(pairs, 'unwrap', geocode=True)
        vert_disp_mm = sbas.vertical_displacement_mm(unwraps_ll)

        Calculate vertical displacement for detrended unwrapped phase grids in geographic coordinates:
        vert_disp_mm = sbas.open_grids(pairs, 'detrend', geocode=True, func=sbas.vertical_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.
        """
        import numpy as np
    
        assert self.is_geo(unwraps), 'ERROR: unwrapped phase defined in radar coordinates'
        
        los_disp = self.los_displacement_mm(unwraps)
        incidence_ll = self.incidence_angle()
        return los_disp/np.cos(incidence_ll)

    def eastwest_displacement_mm(self, unwraps):
        """
        Compute East-West displacement in millimeters.

        Parameters
        ----------
        unwraps : xarray.DataArray or xarray.Dataset
            Unwrapped phase grid(s) in geographic coordinates.

        Returns
        -------
        xarray.DataArray or xarray.Dataset
            East-West displacement grid(s) in millimeters.

        Examples
        --------
        Calculate East-West displacement for unwrapped phase grids in geographic coordinates:
        unwraps_ll = sbas.open_grids(pairs, 'unwrap', geocode=True)
        ew_disp_mm = sbas.eastwest_displacement_mm(unwraps_ll)

        Calculate East-West displacement for detrended unwrapped phase grids in geographic coordinates:
        ew_disp_mm = sbas.open_grids(pairs, 'detrend', geocode=True, func=sbas.eastwest_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.
        """
        import numpy as np
    
        # this displacement is not symmetrical for the orbits due to scene geometries
        orbit = self.df.orbit.unique()[0]
        sign = 1 if orbit == 'D' else -1
        los_disp = self.los_displacement_mm(unwraps)
        incidence_ll = self.incidence_angle()
        return sign * los_disp/np.sin(incidence_ll)

    def sat_look(self, interactive=False):
        import dask
        import xarray as xr
        import numpy as np
        import os
        import sys

        # ..., look_E, look_N, look_U
        satlook_map = {0: 'look_E', 1: 'look_N', 2: 'look_U'}

        def SAT_look(z, lat, lon):
            coords = np.column_stack([lon.ravel(), lat.ravel(), z.ravel()])
            # look_E look_N look_U
            look = self.PRM().SAT_look(coords, binary=True)\
                                     .astype(np.float32)\
                                     .reshape(z.shape[0], z.shape[1], 6)[...,3:]
            return look

        ################################################################################
        # define valid area checking every 10th pixel per the both dimensions
        ################################################################################
        # reference grid
        grid_ll = self.get_intf_ra2ll()
        # do not use coordinate names lat,lon because the output grid saved as (lon,lon) in this case...
        dem = self.get_dem().interp_like(grid_ll).rename({'lat': 'yy', 'lon': 'xx'})
        # prepare lazy coordinate grids
        lat = xr.DataArray(dem.yy.astype(np.float32).chunk(-1))
        lon = xr.DataArray(dem.xx.astype(np.float32).chunk(-1))
        lats, lons = xr.broadcast(lat, lon)
        # unify chunks
        lats = lats.chunk(dem.chunks)
        lons = lons.chunk(dem.chunks)

        # xarray wrapper for the valid area only
        enu = xr.apply_ufunc(
            SAT_look,
            dem,
            lats,
            lons,
            dask='parallelized',
            vectorize=False,
            output_dtypes=[np.float32],
            output_core_dims=[['enu']],
            dask_gufunc_kwargs={'output_sizes': {'enu': 3}}
        )

        # transform to separate variables
        keys_vars = {val: enu[...,key] for (key, val) in satlook_map.items()}
        sat_look = xr.Dataset(keys_vars)

        if interactive:
            return sat_look

        # save to NetCDF file
        filename = self.get_filenames(None, None, 'sat_look')
        if os.path.exists(filename):
            os.remove(filename)
        encoding = {val: self.compression() for (key, val) in satlook_map.items()}
        handler = sat_look.to_netcdf(filename,
                                        encoding=encoding,
                                        engine=self.engine,
                                        compute=False)
        return handler


    def sat_look_parallel(self, interactive=False):
        """
        Build and save satellite look vectors in geographic coordinates.

        Parameters
        ----------
        n_jobs : int, optional
            Number of parallel processing jobs. n_jobs=-1 means using all available processor cores.
        interactive : bool, optional
            If True, returns the delayed computation object for further processing (default is False).

        Returns
        -------
        None or dask.delayed object
            If interactive is False, the function executes the computation and returns None. If interactive is True,
            it returns the delayed computation object for further processing.

        Examples
        --------
        Build and save satellite look vectors:
        sbas.sat_look_parallel()

        Returns
        -------
        None or dask.delayed object
            If interactive is False, the function executes the computation and returns None. If interactive is True,
            it returns the delayed computation object for further processing.

        Notes
        -----
        This function builds and saves satellite look vectors in geographic coordinates. It leverages parallel processing
        using Dask. If interactive is True, it returns the delayed computation object for further processing. Otherwise,
        it executes the computation and returns None.
        """
        import dask

        delayed = self.sat_look(interactive=interactive)

        if not interactive:
            pbar = tqdm_dask(dask.persist(delayed), desc='Satellite Look Vector Computing')
            delayed.compute()
            # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
            import gc; gc.collect()
        else:
            return delayed
