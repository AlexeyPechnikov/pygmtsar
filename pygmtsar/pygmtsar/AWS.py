# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2024, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
class AWS():

    def download_dem(self, geometry, filename=None, product='1s', provider='GLO', n_jobs=4, joblib_backend='loky', skip_exist=True):
        """
        Download Copernicus GLO-30/GLO-90 Digital Elevation Model or NASA SRTM Digital Elevation Model from open AWS storage.

        from pygmtsar import AWS
        dem = AWS().download_dem(AOI)
        dem.plot.imshow()

        AWS().download_dem(S1.scan_slc(DATADIR), 'dem.nc')
        """
        print('NOTE: AWS module is removed. Download DEM using Tiles().download_dem(AOI, filename=...).')
