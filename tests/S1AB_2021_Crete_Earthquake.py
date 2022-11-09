#!/usr/bin/env python3
# S1A_IW_SLC__1SDV_20210929T162245_20210929T162312_039899_04B89A_651D.zip
# S1B_IW_SLC__1SDV_20210923T162204_20210923T162231_028828_0370C4_681D.zip

from pygmtsar import SBAS
import numpy as np
import xarray as xr
from dask.distributed import Client
import matplotlib.pyplot as plt

#Define Parameters
WORKDIR      = 'raw'
DATADIR      = 'data'
CORRLIMIT    = 0.10
DEFOMAX      = 0

if __name__ == '__main__':
    # Run Local Dask Cluster
    client = Client()
    print (client)
    # Init SBAS
    sbas = SBAS(DATADIR, basedir=WORKDIR)
    print (sbas.to_dataframe())
    # Download orbits
    sbas.download_orbits()
    print (sbas.to_dataframe())
    # Reframe scenes
    sbas.set_pins([25.3, 35.0, None, None])
    sbas.reframe_parallel()
    print (sbas.get_pins())
    sbas.to_dataframe()
    # Download DEM
    #sbas.download_dem(backend='GMT', product='SRTM3')
    sbas.download_dem(product='SRTM1')
    print (sbas.get_dem())
    # Align a Stack of Images
    sbas.stack_parallel()
    # DEM in Radar Coordinates
    sbas.topo_ra_parallel()
    print (sbas.get_topo_ra())
    #Interferograms
    pairs = [sbas.to_dataframe().index]
    #decimator = lambda dataarray: dataarray.coarsen({'y': 4, 'x': 4}, boundary='trim').mean()
    # use the DEM default resolution 60 m
    decimator = sbas.pixel_decimator(resolution_meters=60, debug=True)
    sbas.intf_parallel(pairs, func=decimator)
    # geocode matrices
    sbas.geocode_parallel(pairs)
    print (pairs)
    print (sbas.open_grids(pairs, 'phasefilt'))
    print (sbas.open_grids(pairs, 'phasefilt', geocode=True))
    print (sbas.open_grids(pairs, 'corr'))
    # Download landmask
    sbas.download_landmask()
    print (sbas.get_landmask())
    landmask_ll = sbas.get_landmask()
    print (sbas.open_grids(pairs, 'phasefilt', mask=landmask_ll, geocode=True)[0])
    # Unwrapping
    landmask_ra = sbas.get_landmask(inverse_geocode=True)
    cleaner = lambda corr, unwrap: xr.where(corr>=CORRLIMIT, unwrap, np.nan)
    sbas.unwrap_parallel(pairs, threshold=CORRLIMIT, mask=landmask_ra, func=cleaner)
    print (sbas.open_grids(pairs, 'unwrap', geocode=True)[0])
    print (sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm, geocode=True)[0])
    # Incidence angle
    sbas.sat_look_parallel()
    print (sbas.incidence_angle())
    # Displacements
    print (sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm, geocode=True)[0])
    print (sbas.open_grids(pairs, 'unwrap', func=[sbas.vertical_displacement_mm, sbas.nearest_grid], mask=sbas.intf_ra2ll(landmask_ra), geocode=True)[0])
    print (sbas.open_grids(pairs, 'unwrap', func=[sbas.eastwest_displacement_mm, sbas.nearest_grid], mask=sbas.intf_ra2ll(landmask_ra), geocode=True)[0])
