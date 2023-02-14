#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_geocode import SBAS_geocode
from .tqdm_dask import tqdm_dask

class SBAS_incidence(SBAS_geocode):

    def get_sat_look(self):
        import xarray as xr
        import os

        sat_look_file = self.get_filenames(None, None, 'sat_look')
        assert os.path.exists(sat_look_file), 'ERROR: satellite looks grid missed. Build it first using SBAS.sat_look_parallel()'
        sat_look = xr.open_dataset(sat_look_file, engine=self.engine, chunks=self.chunksize).rename({'yy': 'lat', 'xx': 'lon'})

        return sat_look

    #gmt grdmath unwrap_mask.grd $wavel MUL -79.58 MUL = los.grd
    def los_displacement_mm(self, unwraps):
        # constant is negative to make LOS = -1 * range change
        # constant is (1000 mm) / (4 * pi)
        scale = -79.58 * self.PRM().get('radar_wavelength')
        los_disp = scale*unwraps
        return los_disp.rename('los')

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
        encoding = {val: self.compression for (key, val) in satlook_map.items()}
        handler = sat_look.to_netcdf(filename,
                                        encoding=encoding,
                                        engine=self.engine,
                                        compute=False)
        return handler


    def sat_look_parallel(self, interactive=False):
        import dask

        delayed = self.sat_look(interactive=interactive)

        if not interactive:
            tqdm_dask(dask.persist(delayed), desc='Satellite Look Vector Computing')
        else:
            return delayed
