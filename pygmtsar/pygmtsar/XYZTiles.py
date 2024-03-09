# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .datagrid import datagrid
from .tqdm_joblib import tqdm_joblib

class XYZTiles(datagrid, tqdm_joblib):

    http_timeout = 30
    # OSM tiles downloading requires the browser header (otherwise, server returns HTTP 403 code)
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/88.0.4324.150 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
        'Accept-Language': 'en-US,en;q=0.9',
        'Accept-Encoding': 'gzip, deflate, br',
        'Referer': 'https://www.example.com',
        'Connection': 'keep-alive',
        'Cache-Control': 'max-age=0',
        'Upgrade-Insecure-Requests': '1',
        'DNT': '1',
    }

    def download_googlemaps(self, geometry, zoom, filename=None, **kwargs):
        kwargs['url'] = 'https://mt1.google.com/vt/lyrs=r&x={x}&y={y}&z={z}'
        return self.download(geometry, zoom, filename, **kwargs)
        
    def download_googlesatellite(self, geometry, zoom, filename=None, **kwargs):
        kwargs['url'] = 'https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}'
        return self.download(geometry, zoom, filename, **kwargs)

    def download_googlesatellitehybrid(self, geometry, zoom, filename=None, **kwargs):
        kwargs['url'] = 'https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}'
        return self.download(geometry, zoom, filename, **kwargs)

    def download_openstreetmap(self, geometry, zoom, filename=None, **kwargs):
        kwargs['url'] = 'https://tile.openstreetmap.org/{z}/{x}/{y}.png'
        return self.download(geometry, zoom, filename, **kwargs)

    def download(self, geometry, zoom, filename=None, url='https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}', n_jobs=8, skip_exist=True, debug=False):
        """
        Downloads map tiles for a specified geometry and zoom level from a given tile map service.

        Parameters
        ----------
        geometry : object
            The area for which to download map tiles, defined as a Shapely geometry, GeoPandas object, or Xarray object.
            This is typically referred to as an Area of Interest (AOI).
        zoom : int
            The zoom level for the map tiles. Higher zoom levels correspond to higher resolution.
        url : str, optional
            The URL template of the tile map service. The placeholders {x}, {y}, {z} should be present in the URL. 
            Default is Google Satellite Hybrid 'https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}'.
        n_jobs : int, optional
            The number of concurrent download jobs. Default is 8.

        Returns
        -------
        Xarray
            An Xarray object containing the RGB raster data for the downloaded map tiles. This object can be used for further analysis and visualization.

        Examples
        --------
        # Download map tiles at zoom level 10
        from pygmtsar import XYZTiles
        gmap = XYZTiles().download(AOI, zoom = 10)
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
        import os

        if filename is not None and os.path.exists(filename) and skip_exist:
            print ('NOTE: Target file exists, return it. Use "skip_exist=False" or omit the filename to allow new downloading.')
            return xr.open_dataarray(filename, engine=self.netcdf_engine, chunks=self.chunksize)

        def deg2num(lat_deg, lon_deg, zoom):
            lat_rad = math.radians(lat_deg)
            n = 2.0 ** zoom
            x = int((lon_deg + 180.0) / 360.0 * n)
            y = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
            return (x, y)

        def download_tile(url, x, y, z, debug=False):
            #url = f'https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}'
            url_tile = url.format(x=x, y=y, z=z)
            if debug:
                print('DEBUG: XYZTiles: url', url_tile)
            response = requests.get(url_tile, headers=self.headers, timeout=self.http_timeout)
            if response.status_code == 200:
                # Check if the content type is an image
                if 'image' in response.headers.get('Content-Type', ''):
                    image = iio.imread(io.BytesIO(response.content))
                    return image
                else:
                    raise ValueError(f'Expected an image response, got {response.headers.get("Content-Type")}')
            else:
                raise ValueError(f'Request for tile {url_tile} failed with status {response.status_code}')

        def num2deg(xtile, ytile, zoom):
            n = 2.0 ** zoom
            lon_deg = xtile / n * 360.0 - 180.0
            lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
            lat_deg = math.degrees(lat_rad)
            return (lat_deg, lon_deg)

#         # latitudes inverted
#         #lat_start, lat_end = 34.42, 34.36
#         #lon_start, lon_end = -118.5, -118.4
#         if isinstance(geometry, (xr.DataArray, xr.Dataset)):
#             lon_start = geometry.lon.min()
#             lat_start = geometry.lat.min()
#             lon_end   = geometry.lon.max()
#             lat_end   = geometry.lat.max()
#         else:
#             lon_start, lat_start, lon_end, lat_end = geometry.bounds
#         lat_start, lat_end = lat_end, lat_start

        bounds = self.get_bounds(geometry)
        lon_start, lat_start, lon_end, lat_end = bounds
        # latitudes inverted
        lat_start, lat_end = lat_end, lat_start

        # it produces 4 tiles for cases like (39.5, 39.5, 40.0, 40.0)
        #left, right = int(np.floor(lon_start)), int(np.floor(lon_end))
        #bottom, top = int(np.floor(lat_start)), int(np.floor(lat_end))
        # enhancement to produce a single tile for cases like (39.5, 39.5, 40.0, 40.0)
        left = np.floor(min(lon_start, lon_end))
        right = np.ceil(max(lon_start, lon_end)) - 1
        bottom = np.floor(min(lat_start, lat_end))
        top = np.ceil(max(lat_start, lat_end)) - 1
        left, right = int(left), int(right)
        bottom, top = int(bottom), int(top)
        #print ('left, right', left, right, 'bottom, top', bottom, top)
        
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
        def job_tile(x, y, debug=False):
            tile_image = download_tile(url, x, y, zoom, debug)
            tile_array = np.array(tile_image)
            lat_upper_left, lon_upper_left = num2deg(x, y, zoom)
            lat_lower_right, lon_lower_right = num2deg(x + 1, y + 1, zoom)
            # Exclude the endpoint to prevent overlap with the adjacent tile
            latitudes = np.linspace(lat_upper_left, lat_lower_right, tile_array.shape[0], endpoint=False)
            longitudes = np.linspace(lon_upper_left, lon_lower_right, tile_array.shape[1], endpoint=False)
            return xr.DataArray(tile_array, dims=('lat', 'lon', 'band'), 
                                coords={'lat': latitudes, 'lon': longitudes}).rename('colors')

        if n_jobs is None:
            # do not use joblib parallel processing
            tile_xarrays = []
            with self.tqdm_joblib(tqdm(desc=f'XYZ Tiles Downloading', total=(x_end-x_start+1)*(y_end-y_start+1))) as pbar:
                for x in range(x_start, x_end + 1):
                    for y in range(y_start, y_end + 1):
                        tile = job_tile(x, y, debug=debug)
                        tile_xarrays.append(tile)     
                        pbar.update(1)       
        else:
            with self.tqdm_joblib(tqdm(desc='XYZ Tiles Downloading', total=(x_end-x_start+1)*(y_end-y_start+1))) as progress_bar:
                tile_xarrays = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(job_tile)(x, y)\
                                    for x in range(x_start, x_end + 1) for y in range(y_start, y_end + 1))

        da = xr.combine_by_coords(tile_xarrays)['colors'].transpose('band', 'lat', 'lon')
        # fix for inverted latitudes
        if da.lat.diff('lat')[0].item() < 0:
            da = da.reindex(lat=da.lat[::-1])
        # crop geometry extent
        da = da.sel(lat=slice(bounds[1], bounds[3]), lon=slice(bounds[0], bounds[2]))
        
        if filename is not None:
            if os.path.exists(filename):
                os.remove(filename)
            encoding = {'colors': self._compression(da.shape)}
            da.to_netcdf(filename, encoding=encoding, engine=self.netcdf_engine)
        return da
