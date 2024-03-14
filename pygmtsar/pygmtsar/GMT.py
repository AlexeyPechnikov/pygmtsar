# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2024, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
class GMT():
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
        GMT().download_dem(S1.scan_slc(DATADIR))
        
        Notes
        --------
        https://docs.generic-mapping-tools.org/6.0/datasets/earth_relief.html
        """
        print('NOTE: GMT module is removed. Download DEM using Tiles().download_dem(AOI, filename=...).')

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
        landmask = GMT().download_landmask(S1.scan_slc(DATADIR))

        Notes
        -----
        This method downloads the landmask using GMT's local data or server.
        """
        print('NOTE: GMT module is removed. Download landmask using Tiles().download_landmask(AOI, filename=...).')
