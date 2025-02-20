# ----------------------------------------------------------------------------
# InSAR.dev
# 
# This file is part of the InSAR.dev project: https://InSAR.dev
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_sbas import Stack_sbas
from insardev_toolkit import tqdm_dask

class Stack_incidence(Stack_sbas):

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
        los_projection_mm = stack.los_projection(tidal)
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
        stack.los_projection([22.67, 13.36, 0])
        # Expected output:
        # NOTE: estimation using central point satellite look vector
        # array([-15.57419278])
        """
        import xarray as xr
        import numpy as np

        sat_look = self.get_satellite_look_vector()

        if isinstance(data, xr.Dataset):
            los = xr.dot(xr.concat([sat_look.look_E, sat_look.look_N, sat_look.look_U], dim='dim'),
                   xr.concat([data.dx, data.dy, data.dz], dim='dim'),
                  dims=['dim'])
            return los.transpose('date',...)
        elif isinstance(data, (list, tuple)):
            print ('NOTE: estimation using central point satellite look vector')
            look = sat_look.isel(y=sat_look.y.size//2, x=sat_look.x.size//2)
            data = np.column_stack(data)
            return np.dot(data, [look.look_E, look.look_N, look.look_U])

    def get_satellite_look_vector(self):
        """
        Return satellite look vectors in geographic coordinates as Xarray Dataset.

        Returns
        -------
        xarray.Dataset
            The satellite look vectors in geographic coordinates.

        Examples
        --------
        Get satellite look vectors:
        sat_look_ll = stack.get_satellite_look_vector()

        Notes
        -----
        This function returns the satellite look vectors in geographic coordinates as Xarray Dataset. The satellite look vectors
        should be computed and saved prior to calling this function using the `sat_look` method.
        """
        return self.open_cube('sat_look')

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
        unwraps_ra = stack.open_grids(pairs, 'unwrap')
        los_disp_ra = stack.los_displacement_mm(unwraps_ra)
        # or the same code in one line
        los_disp_ra = stack.open_grids(pairs, 'unwrap', func=stack.los_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.

        Calculate LOS displacement for detrended unwrapped phase grids in geographic coordinates:
        detrend_ll = stack.open_grids(pairs, 'detrend', geocode=True)
        los_disp_ll = stack.los_displacement_mm(detrend_ll)
        # or the same code in one line
        los_disp_ll = stack.open_grids(pairs, 'detrend', geocode=True, func=stack.los_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.
        """
        import xarray as xr
        import numpy as np

        # constant is negative to make LOS = -1 * range change
        # constant is (1000 mm) / (4 * pi)
        scale = -79.58 * self.PRM().get('radar_wavelength')

        if isinstance(data, (list, tuple)):
            return scale*np.asarray(data)
        elif isinstance(data, (xr.DataArray)):
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
        inc_angle_ll = stack.incidence_angle()

        Notes
        -----
        This function computes the incidence angle grid in geographic coordinates based on the satellite look vectors.
        The satellite look vectors should be computed and saved prior to calling this function using the `sat_look` method.
        The incidence angle is calculated using the formula:
        incidence_angle = arctan2(sqrt(look_E**2 + look_N**2), look_U)
        """
        import xarray as xr
        import numpy as np

        sat_look = self.get_satellite_look_vector()
        incidence_ll = (np.arctan2(np.sqrt(sat_look.look_E**2 + sat_look.look_N**2), sat_look.look_U) * np.sign(sat_look.look_E))
        return incidence_ll.rename('incidence_angle')

    def plot_incidence_angle(self, data='auto', caption='Incidence Angle in Radar Coordinates, [rad]', cmap='gray', aspect=None, **kwargs):
        import matplotlib.pyplot as plt

        plt.figure()
        if isinstance(data, str) and data == 'auto':
            data = self.incidence_angle()

        data.plot.imshow(cmap=cmap)
        #self.plot_AOI(**kwargs)
        #self.plot_POI(**kwargs)
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        #plt.xlabel('Range')
        #plt.ylabel('Azimuth')
        plt.title(caption)

    def vertical_displacement_mm(self, data):
        """
        Compute vertical displacement in millimeters in radar coordinates.

        Parameters
        ----------
        data : xarray.DataArray or xarray.Dataset
            Unwrapped phase grid(s) in radar coordinates.

        Returns
        -------
        xarray.DataArray
            Vertical displacement grid(s) in millimeters.

        Examples
        --------
        ...
        """
        import numpy as np
        return self.los_displacement_mm(data)/np.cos(self.incidence_angle())

    def eastwest_displacement_mm(self, data):
        """
        Compute East-West displacement in millimeters.

        Parameters
        ----------
        data : xarray.DataArray or xarray.Dataset
            Unwrapped phase grid(s) in geographic coordinates.

        Returns
        -------
        xarray.DataArray or xarray.Dataset
            East-West displacement grid(s) in millimeters.

        Examples
        --------
        Calculate East-West displacement for unwrapped phase grids in geographic coordinates:
        unwraps_ll = stack.open_grids(pairs, 'unwrap', geocode=True)
        ew_disp_mm = stack.eastwest_displacement_mm(unwraps_ll)

        Calculate East-West displacement for detrended unwrapped phase grids in geographic coordinates:
        ew_disp_mm = stack.open_grids(pairs, 'detrend', geocode=True, func=stack.eastwest_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.
        """
        import numpy as np
        return self.los_displacement_mm(data)/np.sin(self.incidence_angle())

    def elevation_m(self, data, baseline=1):
        """
        Computes the elevation in meters from unwrapped phase grids in radar coordinates.

        Parameters
        ----------
        data : xarray.DataArray or xarray.Dataset
            The unwrapped phase grid(s) in radar coordinates.
        
        baseline : numeric, optional
            The perpendicular baseline in meters, by default 1.

        Returns
        -------
        xarray.DataArray
            The elevation grid(s) in meters.

        Examples
        --------
        # Example usage
        elevation_data = sbas.elevation_m(phase_data, baseline=300)
        """
        import xarray as xr
        import numpy as np

        # expected accuracy about 0.01%
        #wavelength, slant_range = self.PRM().get('radar_wavelength','SC_height')
        wavelength, slant_range_start,slant_range_end = self.PRM().get('radar_wavelength', 'SC_height_start', 'SC_height_end')

        incidence_angle = self.incidence_angle()
        slant_range = xr.DataArray(np.linspace(slant_range_start,slant_range_end, incidence_angle.shape[1]),
                                       coords=[incidence_angle.coords['x']])

        if 'stack' in data.dims and 'y' in data.coords and 'x' in data.coords:
            incidence = incidence_angle.interp(y=data.y, x=data.x, method='linear')
            slant = slant_range.interp(x=data.x, method='linear')
        else:
            incidence = incidence_angle.reindex_like(data, method='nearest')
            slant = slant_range.reindex_like(data.x, method='nearest')

        # sign corresponding to baseline and phase signs
        return -(wavelength*data*slant*np.cos(incidence)/(4*np.pi*baseline)).rename('ele')
