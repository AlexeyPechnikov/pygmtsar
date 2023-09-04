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

    def los_projection(self, data):
        """
        Calculate LOS projection for vector defined by its dx, dy, dz components.

        Parameters
        ----------
        data : xarray dataset
            The input data containing the displacement components dx, dy, dz.

        Returns
        -------
        float, numpy.ndarray, pandas.DataFrame
            The LOS projection. Type of return depends on the input type.

        Examples
        -------
        Calculate tidal LOS projection measured in meter [m]:
        los_projection_mm = sbas.los_projection(tidal)
        # Expected input
        # xarray.Dataset
        # Dimensions:
        # date: 31 y: 1 x: 1
        # Data variables:
        # dx (date, y, x) float64 -0.06692 -0.03357 ... 0.005664
        # dy (date, y, x) float64 -0.004765 0.01228 ... -0.04304
        # dz (date, y, x) float64 0.0162 -0.0999 ... 0.005759
        # ...        
        # Expected output:
        # xarray.DataArray date: 31 y: 1 x: 1
        # array([ 0.05532877, -0.05658128, -0.11400223, -0.06658935, -0.0071757 ,
        #    -0.02071992, -0.07211125, -0.12153598, -0.09518547, -0.10037747,
        #    -0.0914933 , -0.12743347, -0.11006747, -0.0643307 , -0.04372583,
        #    -0.07117568, -0.13215618, -0.10467723, -0.01379629,  0.03088265,
        #     0.02786578, -0.01465195, -0.12157386, -0.11801581, -0.001239  ,
        #     0.11614589,  0.07466661, -0.05334002, -0.10686331, -0.06112201,
        #     0.00554765])
        # ...
    
        Calculate plate velocity LOS projection in millimeter [mm]:
        sbas.los_projection([22.67, 13.36, 0])
        # Expected output:
        # NOTE: estimation using central point satellite look vector
        # array([-15.57419278])
        """
        import xarray as xr
        import numpy as np

        sat_look = self.get_sat_look()

        if isinstance(data, xr.Dataset):
            look = sat_look.interp_like(data, method='linear', assume_sorted=True)
            los = xr.dot(xr.concat([look.look_E, look.look_N, look.look_U], dim='dim'),
                   xr.concat([data.dx, data.dy, data.dz], dim='dim'),
                  dims=['dim'])
            return los.transpose('date',...)
        elif isinstance(data, (list, tuple)):
            print ('NOTE: estimation using central point satellite look vector')
            look = sat_look.isel(y=sat_look.y.size//2, x=sat_look.x.size//2)
            data = np.column_stack(data)
            return np.dot(data, [look.look_E, look.look_N, look.look_U])

#     def los_projection(self, data):
#         """
#         Calculate LOS projection for vector defined by its dx, dy, dz components.
# 
#         Parameters
#         ----------
#         data : list, tuple, numpy.ndarray, pandas.DataFrame
#             The input data containing the displacement components dx, dy, dz.
#         scale : float, optional
#             Scale factor to convert input displacements in meter to output LOS displacement in millimeter.
# 
#         Returns
#         -------
#         float, numpy.ndarray, pandas.DataFrame
#             The LOS projection. Type of return depends on the input type.
# 
#         Note
#         ----
#         The function is not optimized for delayed execution.
# 
#         Examples
#         -------
#         Calculate tidal LOS projection:
#         los_projection_mm = 1000*sbas.los_projection(tidal)
#         # Expected input
#         #        lon       lat       dx         dy          dz
#         # date						
#         # 2022-06-16  13400.758  47401.431  -66.917528  -4.765059  16.200381
#         # ...        
#         # Expected output:
#         #        lon       lat       dx         dy          dz        los
#         # date						
#         # 2022-06-16  13400.758  47401.431  -66.917528  -4.765059  16.200381  55.340305
#         # ...
# 
#         Using list or tuple as input:
#         los_projection_mm = 1000*sbas.los_projection([tidal.dx, tidal.dy, tidal.dz], lon, lat)
#         # Expected output:
#         # [55.34030452, -56.55791618, ...]
# 
#         Using numpy.ndarray as input:
#         los_projection_mm = 1000*sbas.los_projection(np.column_stack([tidal.dx, tidal.dy, tidal.dz]))
#         # Expected output (with central point satellite look vector estimation):
#         # [54.72536278, -57.87347137, ...]
# 
#         Note: When lat and lon are not provided, the function will estimate using a central point satellite look vector.
# 
#         """
#         import xarray as xr
#         import pandas as pd
#         import numpy as np
# 
#         sat_look = self.get_sat_look()
# 
#         if isinstance(data, (list, tuple)):
#             data = np.column_stack(data)
# 
#         if isinstance(data, np.ndarray):
#             #if y is not None and x is not None:
#             #    look = sat_look.sel(y=y, x=x, method='nearest')
#             #else:
#             print ('NOTE: estimation using central point satellite look vector')
#             look = sat_look.isel(y=sat_look.y.size//2, x=sat_look.x.size//2)
#             # only for input scalars
#             #return data[0] * sat_look.look_E.values + data[1] * sat_look.look_N.values + data[2] * sat_look.look_U.values
#             return np.dot(data, [look.look_E, look.look_N, look.look_U])
#         elif isinstance(data, pd.DataFrame):
#             # TODO: allow to process multiple coordinates
#             if 'y' in data.columns and 'x' in data.columns:
#                 x = data.loc[data.index[0], 'x']
#                 y = data.loc[data.index[0], 'y']
#                 look = sat_look.sel(y=y, x=x, method='nearest')
#             else:
#                 print ('NOTE: estimation using central point satellite look vector')
#                 look = sat_look.isel(y=sat_look.y.size//2, x=sat_look.x.size//2)
#         elif isinstance(data, xr.Dataset):
#             # TODO: allow to process multiple coordinates
#             if 'y' in data and 'x' in data:
#                 x = data.x.values
#                 y = data.y.values
#                 #print ('y, x', y, x)
#                 look = sat_look.sel(y=y, x=x, method='nearest')
#             else:
#                 print ('NOTE: estimation using central point satellite look vector')
#                 look = sat_look.isel(y=sat_look.y.size//2, x=sat_look.x.size//2)
#         los = np.dot(np.column_stack([data.dx, data.dy, data.dz]), [look.look_E, look.look_N, look.look_U])
#         return los

#     def los_projection(self, data, scale=1):
#         """
#         Calculate LOS projection for vector defined by its dx, dy, dz components.
# 
#         Parameters
#         ----------
#         data : list, tuple, numpy.ndarray, pandas.DataFrame
#             The input data containing the displacement components dx, dy, dz.
#         scale : float, optional
#             Scale factor to convert input displacements in meter to output LOS displacement in millimeter.
# 
#         Returns
#         -------
#         float, numpy.ndarray, pandas.DataFrame
#             The LOS projection. Type of return depends on the input type.
# 
#         Note
#         ----
#         The function is not optimized for delayed execution.
# 
#         Examples
#         -------
#         Calculate tidal LOS projection:
#         los_projection_mm = 1000*sbas.los_projection(tidal)
#         # Expected input
#         #        lon       lat       dx         dy          dz
#         # date						
#         # 2022-06-16  13400.758  47401.431  -66.917528  -4.765059  16.200381
#         # ...        
#         # Expected output:
#         #        lon       lat       dx         dy          dz        los
#         # date						
#         # 2022-06-16  13400.758  47401.431  -66.917528  -4.765059  16.200381  55.340305
#         # ...
# 
#         Using list or tuple as input:
#         los_projection_mm = 1000*sbas.los_projection([tidal.dx, tidal.dy, tidal.dz], lon, lat)
#         # Expected output:
#         # [55.34030452, -56.55791618, ...]
# 
#         Using numpy.ndarray as input:
#         los_projection_mm = 1000*sbas.los_projection(np.column_stack([tidal.dx, tidal.dy, tidal.dz]))
#         # Expected output (with central point satellite look vector estimation):
#         # [54.72536278, -57.87347137, ...]
# 
#         Note: When lat and lon are not provided, the function will estimate using a central point satellite look vector.
# 
#         """
#         import pandas as pd
#         import numpy as np
# 
#         sat_look = self.get_sat_look()
# 
#         if isinstance(data, (list, tuple)):
#             data = np.column_stack(data)
# 
#         if isinstance(data, np.ndarray):
#             #if y is not None and x is not None:
#             #    look = sat_look.sel(y=y, x=x, method='nearest')
#             #else:
#             print ('NOTE: estimation using central point satellite look vector')
#             look = sat_look.isel(y=sat_look.y.size//2, x=sat_look.x.size//2)
#             # only for input scalars
#             #return data[0] * sat_look.look_E.values + data[1] * sat_look.look_N.values + data[2] * sat_look.look_U.values
#             return np.dot(data, [look.look_E, look.look_N, look.look_U])
# 
#         #elif isinstance(data, pd.DataFrame):
#         # TODO: allow to process multiple coordinates
#         if 'y' in data.columns and 'x' in data.columns:
#             x = data.loc[data.index[0], 'x']
#             y = data.loc[data.index[0], 'y']
#             look = sat_look.sel(y=y, x=x, method='nearest')
#         else:
#             print ('NOTE: estimation using central point satellite look vector')
#             look = sat_look.isel(y=sat_look.y.size//2, x=sat_look.x.size//2)
#         los = np.dot(np.column_stack([data.dx, data.dy, data.dz]), [look.look_E, look.look_N, look.look_U])
#         return data.assign(los=scale*los)

    def get_sat_look(self, chunksize=None):
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
        subswath = self.get_subswath()    
        return self.open_grid('sat_look', subswath=subswath, chunksize=chunksize)

    #gmt grdmath unwrap_mask.grd $wavel MUL -79.58 MUL = los.grd
    def los_displacement_mm(self, data):
        """
        Compute line-of-sight (LOS) displacement in millimeters.

        Parameters
        ----------
        data : xarray.DataArray or constant, list, tuple, Numpy array, Pandas Series
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
        import xarray as xr
        import numpy as np

        # constant is negative to make LOS = -1 * range change
        # constant is (1000 mm) / (4 * pi)
        scale = -79.58 * self.PRM().get('radar_wavelength')

        if isinstance(data, (list, tuple)):
            print ('X')
            return scale*np.asarray(data)
        elif isinstance(data, (xr.DataArray)):
            print ('Y')
            return (scale*data).rename('los')
        else:
            return scale*data

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

    def sat_look_parallel(self, chunksize=None, interactive=False):
        #import dask
        import xarray as xr
        import numpy as np

        subswath = self.get_subswath()

        # ..., look_E, look_N, look_U
        satlook_map = {0: 'look_E', 1: 'look_N', 2: 'look_U'}

        def SAT_look(z, lat, lon):
            coords = np.column_stack([lon.ravel(), lat.ravel(), z.ravel()])
            # look_E look_N look_U
            look = self.PRM().SAT_look(coords, binary=True)\
                                     .astype(np.float32)\
                                     .reshape(z.shape[0], z.shape[1], 6)[...,3:]
            return look

        if chunksize is None:
            chunksize = self.chunksize

        # reference grid
        trans_inv = self.get_trans_inv(subswath, chunksize=chunksize)[['lt', 'll', 'ele']]

        # xarray wrapper for the valid area only
        enu = xr.apply_ufunc(
            SAT_look,
            trans_inv.ele,
            trans_inv.lt,
            trans_inv.ll,
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

        return self.save_grid(sat_look, 'sat_look', subswath,
                              f'Satellite Look Vector Computing sw{subswath}', chunksize)
