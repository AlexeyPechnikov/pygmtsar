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
    # https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N38_00_E038_00_DEM/Copernicus_DSM_COG_10_N38_00_E038_00_DEM.tif
    base_url_glo = 'https://copernicus-dem-{resolution}m.s3.amazonaws.com'
    path_id_glo = 'Copernicus_DSM_COG_{product1}0_{SN2}_00_{WE3}_00_DEM'
    # Copernicus_DSM_COG_10_S16_00_W043_00_DEM
    tile_id_glo = path_id_glo + '.tif'

    #aws s3 ls --no-sign-request s3://elevation-tiles-prod/skadi/
    # https://s3.amazonaws.com/elevation-tiles-prod/skadi/N20/N20E000.hgt.gz
    base_url_srtm = 'https://s3.amazonaws.com/elevation-tiles-prod/skadi'
    path_id_srtm = '{SN2}'
    tile_id_srtm = '{SN2}{WE3}.hgt.gz'

    def _download_tile_glo(self, product, lon, lat):
        """
        Download Copernicus GLO-30/GLO-90 Digital Elevation Model tiles from open AWS storage.
        
        Previously, it was in-memory processing but it seems not stable on slow internet connection:
        import io
        with io.BytesIO(response.content) as f:
            #tile = xr.open_dataarray(f, engine='rasterio')...
            tile = rio.open_rasterio(f, chunks=self.chunksize)...
        """
        import rioxarray as rio
        import requests
        import os
        import tempfile

        product1 = int(product[0])
        SN2 = f'{"S" if lat<0 else "N"}{abs(lat):02}'
        WE3 = f'{"W" if lon<0 else "E"}{abs(lon):03}'
        base_url = self.base_url_glo.format(resolution=30 if product=='1s' else 90)
        path_id = self.path_id_glo.format(product1=product1, SN2=SN2, WE3=WE3)
        tile_id = self.tile_id_glo.format(product1=product1, SN2=SN2, WE3=WE3)
        tile_url = f'{base_url}/{path_id}/{tile_id}'
        tile_filename = os.path.join(tempfile.gettempdir(), tile_id)
        print ('tile_url', tile_url)
        print ('tile_filename', tile_filename)
        try:
            with requests.get(tile_url, stream=True) as response:
                response.raise_for_status()
                with open(tile_filename, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=256*1024):
                        f.write(chunk)
            with rio.open_rasterio(tile_filename) as raster:
                tile = raster\
                    .squeeze(drop=True)\
                    .rename({'y': 'lat', 'x': 'lon'})\
                    .drop_vars('spatial_ref')\
                    .load()
            if tile.lat.diff('lat')[0].item() < 0:
                tile = tile.reindex(lat=tile.lat[::-1])
        except requests.exceptions.RequestException as e:
            # offshore tiles are missed by design
            print(f'Request error for {tile_id}: {e}')
            tile = None
        except Exception as e:
            print(e)
            raise
        finally:
            if os.path.exists(tile_filename):
                os.remove(tile_filename)
        return tile

    def _download_tile_srtm(self, product, lon, lat):
        """
        Download NASA SRTM Digital Elevation Model tiles from open AWS storage.
        """
        import rioxarray as rio
        import requests
        import os
        import tempfile
        import gzip

        SN2 = f'{"S" if lat<0 else "N"}{abs(lat):02}'
        WE3 = f'{"W" if lon<0 else "E"}{abs(lon):03}'
        base_url = self.base_url_srtm
        path_id = self.path_id_srtm.format(SN2=SN2, WE3=WE3)
        tile_id = self.tile_id_srtm.format(SN2=SN2, WE3=WE3)
        tile_url = f'{base_url}/{path_id}/{tile_id}'
        # remove .gz extension
        tile_filename = os.path.join(tempfile.gettempdir(), tile_id[:-3])
        print ('tile_url', tile_url)
        print ('tile_filename', tile_filename)
        try:
            with requests.get(tile_url, stream=True) as response:
                response.raise_for_status()
                # Stream the content under the context of gzip decompression
                with gzip.GzipFile(fileobj=response.raw) as gz:
                    with open(tile_filename, 'wb') as f:
                        while True:
                            chunk = gz.read(256 * 1024)
                            if not chunk:
                                break
                            f.write(chunk)
            with rio.open_rasterio(tile_filename) as raster:
                tile = raster\
                    .squeeze(drop=True)\
                    .rename({'y': 'lat', 'x': 'lon'})\
                    .drop_vars('spatial_ref')\
                    .load()
            if tile.lat.diff('lat')[0].item() < 0:
                tile = tile.reindex(lat=tile.lat[::-1])
        except requests.exceptions.RequestException as e:
            # offshore tiles are missed by design
            print(f'Request error for {tile_id}: {e}')
            tile = None
        except Exception as e:
            print(e)
            raise
        finally:
            if os.path.exists(tile_filename):
                os.remove(tile_filename)
        return tile

    def download_dem(self, geometry, filename=None, product='1s', provider='GLO', n_jobs=4, joblib_backend='loky', skip_exist=True):
        """
        Download Copernicus GLO-30/GLO-90 Digital Elevation Model or NASA SRTM Digital Elevation Model from open AWS storage.

        from pygmtsar import AWS
        dem = AWS().download_dem(AOI)
        dem.plot.imshow()

        AWS().download_dem(S1.scan_slc(DATADIR), 'dem.nc')
        """
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        import os

        assert provider in ['GLO', 'SRTM'], f'ERROR: provider name is invalid: {provider}. Expected names are "GLO", "SRTM".'
        if provider == 'SRTM':
            assert product in ['1s'], f'ERROR: only product="1s" is supported for provider {provider}.'
        if provider == 'GLO':
            assert product in ['1s', '3s'], f'ERROR: product name is invalid: {product} for provider {provider}. Expected names are "1s", "3s".'

        if filename is not None and os.path.exists(filename) and skip_exist:
            print ('NOTE: DEM file exists, ignore the command. Use "skip_exist=False" or omit the filename to allow new downloading.')
            return

        bounds = self.get_bounds(geometry)

        # it produces 4 tiles for cases like (39.5, 39.5, 40.0, 40.0)
        #left, right = int(np.floor(lon_start)), int(np.floor(lon_end))
        #bottom, top = int(np.floor(lat_start)), int(np.floor(lat_end))
        # enhancement to produce a single tile for cases like (39.5, 39.5, 40.0, 40.0)
        lon_start, lat_start, lon_end, lat_end = bounds
        left = np.floor(min(lon_start, lon_end))
        right = np.ceil(max(lon_start, lon_end)) - 1
        bottom = np.floor(min(lat_start, lat_end))
        top = np.ceil(max(lat_start, lat_end)) - 1
        left, right = int(left), int(right)
        bottom, top = int(bottom), int(top)
        #print ('left, right', left, right, 'bottom, top', bottom, top)

        job_tile = self._download_tile_glo if provider == 'GLO' else self._download_tile_srtm
        with self.tqdm_joblib(tqdm(desc=f'AWS {provider} DEM Tiles Downloading', total=(right-left+1)*(top-bottom+1))) as progress_bar:
            tile_xarrays = joblib.Parallel(n_jobs=n_jobs, backend=joblib_backend)(joblib.delayed(job_tile)(product, x, y)\
                                for x in range(left, right + 1) for y in range(bottom, top + 1))

        dem = xr.combine_by_coords([tile for tile in tile_xarrays if tile is not None])
        dem = dem.sel(lat=slice(bounds[1], bounds[3]), lon=slice(bounds[0], bounds[2]))

        if filename is not None:
            if os.path.exists(filename):
                os.remove(filename)
            encoding = {'dem': self._compression(dem.shape)}
            dem.rename('dem').to_netcdf(filename, encoding=encoding, engine=self.netcdf_engine)
        else:
            return dem
