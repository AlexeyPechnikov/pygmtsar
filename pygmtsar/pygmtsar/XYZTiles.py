# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .tqdm_joblib import tqdm_joblib

class XYZTiles(tqdm_joblib):

    def download(self, geometry, zoom, n_jobs=-1, url='https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}'):
        """
        from pygmtsar import XYZTiles
        gmap = XYZTiles().download(AOI.geometry[0], zoom = 17)
        gmap.plot.imshow()
        """
        import xarray as xr
        import requests
        import imageio.v3 as iio
        import math
        import io
        import numpy as np
        from tqdm.auto import tqdm
        import joblib

        def deg2num(lat_deg, lon_deg, zoom):
            lat_rad = math.radians(lat_deg)
            n = 2.0 ** zoom
            x = int((lon_deg + 180.0) / 360.0 * n)
            y = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
            return (x, y)

        def download_tile(url, x, y, z):
            #url = f'https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}'
            response = requests.get(url.format(x=x, y=y, z=z))
            image = iio.imread(io.BytesIO(response.content))
            return image

        def num2deg(xtile, ytile, zoom):
            n = 2.0 ** zoom
            lon_deg = xtile / n * 360.0 - 180.0
            lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
            lat_deg = math.degrees(lat_rad)
            return (lat_deg, lon_deg)

        # latitudes inverted
        #lat_start, lat_end = 34.42, 34.36
        #lon_start, lon_end = -118.5, -118.4
        if isinstance(geometry, (xr.DataArray, xr.Dataset)):
            lon_start = geometry.lon.min()
            lat_start = geometry.lat.min()
            lon_end   = geometry.lon.max()
            lat_end   = geometry.lat.max()
        else:
            lon_start, lat_start, lon_end, lat_end = geometry.bounds
        lat_start, lat_end = lat_end, lat_start

        # Calculate tile range
        x_start, y_start = deg2num(lat_start, lon_start, zoom)
        x_end, y_end = deg2num(lat_end, lon_end, zoom)

    #     def job_tile(x, y):
    #         tile_image = download_tile(url, x, y, zoom)
    #         tile_array = np.array(tile_image)
    #         lat_upper_left, lon_upper_left = num2deg(x, y, zoom)
    #         lat_lower_right, lon_lower_right = num2deg(x + 1, y + 1, zoom)
    #         latitudes = np.linspace(lat_upper_left, lat_lower_right, tile_array.shape[0])
    #         longitudes = np.linspace(lon_upper_left, lon_lower_right, tile_array.shape[1])
    #         return xr.DataArray(tile_array, dims=('lat', 'lon', 'band'), 
    #                                coords={'lat': latitudes, 'lon': longitudes}).rename('tile')
        def job_tile(x, y):
            tile_image = download_tile(url, x, y, zoom)
            tile_array = np.array(tile_image)
            lat_upper_left, lon_upper_left = num2deg(x, y, zoom)
            lat_lower_right, lon_lower_right = num2deg(x + 1, y + 1, zoom)
            # Exclude the endpoint to prevent overlap with the adjacent tile
            latitudes = np.linspace(lat_upper_left, lat_lower_right, tile_array.shape[0], endpoint=False)
            longitudes = np.linspace(lon_upper_left, lon_lower_right, tile_array.shape[1], endpoint=False)
            return xr.DataArray(tile_array, dims=('lat', 'lon', 'band'), 
                                coords={'lat': latitudes, 'lon': longitudes}).rename('colors')

        with self.tqdm_joblib(tqdm(desc='XYZ Tile Downloading', total=(x_end-x_start+1)*(y_end-y_start+1))) as progress_bar:
            tile_xarrays = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(job_tile)(x, y)\
                                for x in range(x_start, x_end + 1) for y in range(y_start, y_end + 1))

        return xr.combine_by_coords(tile_xarrays)['colors'].transpose('band', 'lat', 'lon')
