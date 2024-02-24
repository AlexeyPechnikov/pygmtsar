# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2024, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .datagrid import datagrid
from .tqdm_joblib import tqdm_joblib

class AWS(datagrid, tqdm_joblib):
    # https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com
    base_url = 'https://copernicus-dem-{resolution}m.s3.amazonaws.com'
    # Copernicus_DSM_COG_10_S16_00_W043_00_DEM
    tile_id = 'Copernicus_DSM_COG_{product1}0_{SN2}_00_{WE3}_00_DEM'

    def download_dem(self, geometry, filename=None, n_jobs=-1, product='1s', skip_exist=True):
        """
        Download Copernicus GLO-30/GLO-90 Digital Elevation Model from open AWS storage.

        from pygmtsar import AWS
        dem = AWS().download_dem(AOI.geometry[0])
        dem.plot.imshow()
        """
        import xarray as xr
        import rioxarray
        import geopandas as gpd
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        import requests
        import io
        import os

        assert product in ['1s', '3s'], f'ERROR: product name is invalid: {product}. Expected names are "1s", "3s"'

        if filename is not None and os.path.exists(filename) and skip_exist:
            print ('NOTE: DEM file exists, ignore the command. Use "skip_exist=False" or omit the filename to allow new downloading')
            return

        lon_start, lat_start, lon_end, lat_end = self.get_bounds(geometry)

        def job_tile(product, lon, lat):
            base_url = self.base_url.format(resolution=30 if product=='1s' else 90)
            tile_id = self.tile_id.format(product1=int(product[0]),
                                          SN2=f'{"S" if lat<0 else "N"}{lat:02}',
                                          WE3=f'{"W" if lon<0 else "E"}{lon:03}')
            url = f'{base_url}/{tile_id}/{tile_id}.tif'
            response = requests.get(url)
            #response.raise_for_status()
            # offshore tiles are missed by design
            if response.status_code != 200:
                return None
            with io.BytesIO(response.content) as f:
                tile = rioxarray.open_rasterio(f, cache=True)\
                    .squeeze(drop=True)\
                    .rename({'y': 'lat', 'x': 'lon'})\
                    .drop_vars('spatial_ref')\
                    .load()
                if tile.lat.diff('lat')[0].item() < 0:
                    tile = tile.reindex(lat=tile.lat[::-1])
            return tile

        left, right = int(np.floor(lon_start)), int(np.floor(lon_end))
        lower, upper = int(np.floor(lat_start)), int(np.floor(lat_end))
        #print ('left, right', left, right, 'lower, upper', lower, upper)
        with self.tqdm_joblib(tqdm(desc='DEM Tile Downloading', total=(right-left+1)*(upper-lower+1))) as progress_bar:
            tile_xarrays = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(job_tile)(product, x, y)\
                                for x in range(left, right + 1) for y in range(lower, upper + 1))
        dem = xr.combine_by_coords([tile for tile in tile_xarrays if tile is not None])

        if filename is not None:
            if os.path.exists(filename):
                os.remove(filename)
            encoding = {'dem': self._compression(dem.shape)}
            dem.rename('dem').to_netcdf(filename, encoding=encoding, engine=self.netcdf_engine)
        else:
            return dem
