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

class GMT(datagrid, tqdm_joblib):

#     def download_dem(self, geometry, filename=None, product='1s', skip_exist=True):
#         """
#         Download and preprocess SRTM digital elevation model (DEM) data using GMT library.
# 
#         Parameters
#         ----------
#         product : str, optional
#             Product type of the DEM data. Available options are '1s' or 'SRTM1' (1 arcsec ~= 30m, default)
#             and '3s' or 'SRTM3' (3 arcsec ~= 90m).
# 
#         Returns
#         -------
#         None or Xarray Dataarray
# 
#         Examples
#         --------
#         Download default STRM1 DEM (~30 meters):
#         GMT().download_dem()
# 
#         Download STRM3 DEM (~90 meters):
#         GMT.download_dem(product='STRM3')
#         
#         Download default STRM DEM to cover the selected area AOI:
#         GMT().download_dem(AOI)
#         
#         Download default STRM DEM to cover all the scenes:
#         GMT().download_dem(stack.get_extent().buffer(0.1))
# 
#         import pygmt
#         # Set the GMT data server limit to N Mb to allow for remote downloads
#         pygmt.config(GMT_DATA_SERVER_LIMIT=1e6)
#         GMT().download_dem(AOI, product='1s')
#         """
#         import numpy as np
#         import pygmt
#         # suppress warnings
#         pygmt.config(GMT_VERBOSE='errors')
#         import os
#         from tqdm.auto import tqdm
# 
#         if product in ['SRTM1', '1s', '01s']:
#             resolution = '01s'
#         elif product in ['SRTM3', '3s', '03s']:
#             resolution = '03s'
#         else:
#             print (f'ERROR: unknown product {product}. Available only SRTM1 ("01s") and SRTM3 ("03s") DEM using GMT servers')
#             return
# 
#         if filename is not None and os.path.exists(filename) and skip_exist:
#             print ('NOTE: DEM file exists, ignore the command. Use "skip_exist=False" or omit the filename to allow new downloading')
#             return
# 
#         lon_start, lat_start, lon_end, lat_end = self.get_bounds(geometry)
#         with tqdm(desc='GMT SRTM DEM Downloading', total=1) as pbar:
#             # download DEM using GMT extent W E S N
#             dem = pygmt.datasets.load_earth_relief(resolution=resolution, region=[lon_start, lon_end, lat_start, lat_end])
#             pbar.update(1)
# 
#         if filename is not None:
#             if os.path.exists(filename):
#                 os.remove(filename)
#             encoding = {'dem': self._compression(dem.shape)}
#             dem.rename('dem').load().to_netcdf(filename, encoding=encoding, engine=self.netcdf_engine)
#         else:
#             return dem

#     def download_landmask(self, geometry, filename=None, product='1s', resolution='f', skip_exist=True):
#         """
#         Download the landmask and save as NetCDF file.
# 
#         Parameters
#         ----------
#         product : str, optional
#                 Available options are '1s' (1 arcsec ~= 30m, default) and '3s' (3 arcsec ~= 90m).
# 
#         Examples
#         --------
#         from pygmtsar import GMT
#         landmask = GMT().download_landmask(stack.get_dem())
# 
#         Notes
#         -----
#         This method downloads the landmask using GMT's local data or server.
#         """
#         import geopandas as gpd
#         import numpy as np
#         import pygmt
#         # suppress warnings
#         pygmt.config(GMT_VERBOSE='errors')
#         import os
#         #import subprocess
#         from tqdm.auto import tqdm
# 
#         if filename is not None and os.path.exists(filename) and skip_exist:
#             print ('NOTE: landmask file exists, ignore the command. Use "skip_exist=False" or omit the filename to allow new downloading')
#             return
# 
#         if not product in ['1s', '3s']:
#             print (f'ERROR: unknown product {product}. Available only "1s" or "3s" land mask using GMT servers')
#             return
# 
#         lon_start, lat_start, lon_end, lat_end = self.get_bounds(geometry)
#         with tqdm(desc='GMT Landmask Downloading', total=1) as pbar:
#             landmask = pygmt.grdlandmask(resolution=resolution, region=[lon_start, lon_end, lat_start, lat_end], spacing=product, maskvalues='NaN/1')
#             pbar.update(1)
# 
#         if filename is not None:
#             if os.path.exists(filename):
#                 os.remove(filename)
#             encoding = {'landmask': self._compression(landmask.shape)}
#             landmask.rename('landmask').load().to_netcdf(filename, encoding=encoding, engine=self.netcdf_engine)
#         else:
#             return landmask

    def download_dem(self, geometry, filename=None, product='1s', skip_exist=True):
        """
        Download and preprocess SRTM digital elevation model (DEM) data using GMT library.

        Parameters
        ----------
        product : str, optional
            Product type of the DEM data. Available options are '1s' or 'SRTM1' (1 arcsec ~= 30m, default)
            and '3s' or 'SRTM3' (3 arcsec ~= 90m).

        Returns
        -------
        None or Xarray Dataarray

        Examples
        --------
        Download default STRM1 DEM (~30 meters):
        GMT().download_dem()

        Download STRM3 DEM (~90 meters):
        GMT.download_dem(product='STRM3')
    
        Download default STRM DEM to cover the selected area AOI:
        GMT().download_dem(AOI)
    
        Download default STRM DEM to cover all the scenes:
        GMT().download_dem(stack.get_extent().buffer(0.1))
        
        Notes
        --------
        https://docs.generic-mapping-tools.org/6.0/datasets/earth_relief.html
        """
        import os
        import subprocess
        from tqdm.auto import tqdm
        import tempfile
        import warnings
        warnings.filterwarnings('ignore')

        def load_earth_relief(bounds, product, filename):
            if os.path.exists(filename):
                os.remove(filename)
            argv = ['gmt', 'grdcut', f'@earth_relief_{product}', f'-R{bounds[0]}/{bounds[2]}/{bounds[1]}/{bounds[3]}', f'-G{filename}']
            #print ('gmt grdcut argv:', ' '.join(argv))
            p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
            stdout_data, stderr_data = p.communicate()
            if p.returncode != 0:
                print(f'Error executing gmt grdcut: {stderr_data}')
            return stdout_data.strip()

        if product in ['SRTM1', '1s', '01s']:
            resolution = '01s'
        elif product in ['SRTM3', '3s', '03s']:
            resolution = '03s'
        else:
            # use the specified string as is
            pass

        if filename is not None and os.path.exists(filename) and skip_exist:
            print ('NOTE: DEM file exists, ignore the command. Use "skip_exist=False" or omit the filename to allow new downloading')
            return

        grdname = filename
        if grdname is None:
            with tempfile.NamedTemporaryFile() as tmpfile:
                grdname = tmpfile.name + '.grd'

        with tqdm(desc=f'GMT SRTM {product} DEM Downloading', total=1) as pbar:
            load_earth_relief(self.get_bounds(geometry), resolution, grdname)
            pbar.update(1)

        if filename is None:
            da = xr.open_dataarray(grdname, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
            os.remove(grdname)
            return da

    def download_landmask(self, geometry, filename=None, product='1s', resolution='f', skip_exist=True):
        """
        Download the landmask and save as NetCDF file.

        Parameters
        ----------
        product : str, optional
                Available options are '1s' (1 arcsec ~= 30m, default) and '3s' (3 arcsec ~= 90m).

        Examples
        --------
        from pygmtsar import GMT
        landmask = GMT().download_landmask(stack.get_dem())

        Notes
        -----
        This method downloads the landmask using GMT's local data or server.
        """
        import os
        import subprocess
        from tqdm.auto import tqdm
        import tempfile
        import warnings
        warnings.filterwarnings('ignore')

        def grdlandmask(bounds, product, resolution, filename):
            if os.path.exists(filename):
                os.remove(filename)
            argv = ['gmt', 'grdlandmask', f'-R{bounds[0]}/{bounds[2]}/{bounds[1]}/{bounds[3]}', f'-I{product}', f'-D{resolution}', '-NNaN/1', f'-G{filename}']
            #print ('grdlandmask argv:', ' '.join(argv))
            p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
            stdout_data, stderr_data = p.communicate()
            if p.returncode != 0:
                print(f'Error executing grdlandmask: {stderr_data}')
            return stdout_data.strip()

        if filename is not None and os.path.exists(filename) and skip_exist:
            print ('NOTE: landmask file exists, ignore the command. Use "skip_exist=False" or omit the filename to allow new downloading')
            return

        if not product in ['1s', '3s']:
            print (f'ERROR: unknown product {product}. Available only "1s" or "3s" land mask using GMT servers')
            return

        grdname = filename
        if grdname is None:
            with tempfile.NamedTemporaryFile() as tmpfile:
                grdname = tmpfile.name + '.grd'

        with tqdm(desc='GMT Landmask Downloading', total=1) as pbar:
            grdlandmask(self.get_bounds(geometry), product, resolution, grdname)
            pbar.update(1)

        if filename is None:
            da = xr.open_dataarray(grdname, engine=self.netcdf_engine, chunks=self.netcdf_chunksize)
            os.remove(grdname)
            return da
