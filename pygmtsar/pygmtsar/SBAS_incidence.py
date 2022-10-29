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

    def sat_look_parallel(self, n_jobs=-1, interactive=False):
        import numpy as np
        import xarray as xr
        from tqdm.auto import tqdm
        import joblib
        import os

        def SAT_look(ilat, ilon):
            # lazy dask arrays
            lats, lons = xr.broadcast(coordlat[ilat], coordlon[ilon])
            data = grid_ll.sel(lat=xr.DataArray(lats.values.ravel()),
                                 lon=xr.DataArray(lons.values.ravel()),
                                 method='nearest').compute()
            coords = np.column_stack([lons.values.ravel(), lats.values.ravel(), data.values.ravel()])
            # look_E look_N look_U
            look = self.PRM().SAT_look(coords, binary=True)\
                                     .reshape(-1,6)[:,3:].astype(np.float32)
            # prepare output as xarray dataset
            dims = ['lat', 'lon']
            coords = coords={'lat': coordlat[ilat], 'lon':coordlon[ilon]}
            look_E = xr.DataArray(look[:,0].reshape(lats.shape), dims=dims, coords=coords, name='look_E')
            look_N = xr.DataArray(look[:,1].reshape(lats.shape), dims=dims, coords=coords, name='look_N')
            look_U = xr.DataArray(look[:,2].reshape(lats.shape), dims=dims, coords=coords, name='look_U')
            return xr.merge([look_E, look_N, look_U])

        # take elevation values on the interferogram area in geographic coordinates
        grid_ll = self.get_intf_ra2ll()
        dem = self.get_dem()
        grid_ll = dem.interp_like(grid_ll)
    
        lats, lons = grid_ll.data.numblocks
        latchunks, lonchunks = grid_ll.chunks
        coordlat = np.array_split(grid_ll.lat, np.cumsum(latchunks))
        coordlon = np.array_split(grid_ll.lon, np.cumsum(lonchunks))
        with self.tqdm_joblib(tqdm(desc=f'SAT_look Computing', total=lats*lons)) as progress_bar:
            sat_look = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(SAT_look)(ilat, ilon) \
                                           for ilat in range(lats) for ilon in range(lons))
        # concatenate the chunks
        if xr.__version__ == '0.19.0':
            # for Google Colab
            sat_look = xr.merge(sat_look)
        else:
            # for modern xarray versions
            sat_look = xr.combine_by_coords(sat_look)

        # fill NODATA
        sat_look = xr.where(grid_ll != self.noindex, sat_look, np.nan)
    
        if interactive:
            return sat_look

        # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
        #ds.attrs['node_offset'] = 1
        # save to NetCDF file
        sat_look_file = self.get_filenames(None,None,'sat_look')
        # cleanup before creating the new file
        if os.path.exists(sat_look_file):
            os.remove(sat_look_file)
        sat_look.to_netcdf(sat_look_file,
                     encoding={var:self.compression for var in sat_look.data_vars},
                     engine=self.engine)
