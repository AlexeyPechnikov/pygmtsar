#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_geocode import SBAS_geocode

class SBAS_incidence(SBAS_geocode):

    def get_sat_look(self):
        import xarray as xr
        import os

        sat_look_file = self.get_filenames(None, None, 'sat_look')
        assert os.path.exists(sat_look_file), 'ERROR: satellite looks grid missed. Build it first using SBAS.sat_look_parallel()'
        sat_look = xr.open_dataset(sat_look_file, engine=self.engine, chunks=self.chunksize)

        return sat_look

    #gmt grdmath unwrap_mask.grd $wavel MUL -79.58 MUL = los.grd
    def los_displacement_mm(self, unwraps):
        # constant is negative to make LOS = -1 * range change
        # constant is (1000 mm) / (4 * pi)
        scale = -79.58 * self.PRM().get('radar_wavelength')
        los_disp = scale*unwraps
        return los_disp

    def incidence_angle(self):
        import xarray as xr
        import numpy as np

        sat_look = self.get_sat_look()
        incidence_ll = np.arctan2(np.sqrt(sat_look.look_E**2 + sat_look.look_N**2), sat_look.look_U).rename('incidence_angle')
        return incidence_ll

    def vertical_displacement_mm(self, unwraps):
        import numpy as np
    
        assert self.is_geo(unwraps), 'ERROR: unwrapped phase defined in radar coordinates'
        
        los_disp = self.los_displacement_mm(unwraps)
        incidence_ll = self.incidence_angle()
        return los_disp/np.cos(incidence_ll)

    def eastwest_displacement_mm(self, unwraps):
        import numpy as np
    
        # this displacement is not symmetrical for the orbits due to scene geometries
        orbit = self.df.orbit.unique()[0]
        sign = 1 if orbit == 'D' else -1
        los_disp = self.los_displacement_mm(unwraps)
        incidence_ll = self.incidence_angle()
        return sign * los_disp/np.sin(incidence_ll)