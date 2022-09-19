#!/usr/bin/env python3

from pygmtsar import SBAS
import numpy as np
import xarray as xr

#Define Parameters
WORKDIR      = 'raw'
DATADIR      = 'data'
CORRLIMIT    = 0.10
DEFOMAX      = 0

# Init SBAS
sbas = SBAS(DATADIR, basedir=WORKDIR)
print (sbas.to_dataframe())
# Download orbits
sbas.download_orbits()
print (sbas.to_dataframe())
# Download DEM
#sbas.download_dem(backend='GMT')
sbas.download_dem(product='SRTM1')
print (sbas.get_dem())
# Align a Stack of Images
sbas.stack_parallel()
# DEM in Radar Coordinates
sbas.topo_ra_parallel()
print (sbas.get_topo_ra())
#Interferograms
pairs = [sbas.to_dataframe().index]
decimator = lambda dataarray: dataarray.coarsen({'y': 4, 'x': 4}, boundary='trim').mean()
sbas.intf_parallel(pairs, func=decimator)
print (pairs)
print (sbas.open_grids(pairs, 'phasefilt'))
print (sbas.open_grids(pairs, 'phasefilt', geocode=True))
print (sbas.open_grids(pairs, 'corr'))
# Download landmask
sbas.download_landmask()
print (sbas.get_landmask())
# Unwrapping
landmask_ra = sbas.get_landmask(inverse_geocode=True)
cleaner = lambda corr, unwrap: xr.where(corr>=CORRLIMIT, unwrap, np.nan)
conf = sbas.snaphu_config(NTILEROW=2, NTILECOL=2, ROWOVRLP=200, COLOVRLP=200)
sbas.unwrap_parallel(pairs, threshold=CORRLIMIT, mask=landmask_ra, func=cleaner, conf=conf)
print (sbas.open_grids(pairs, 'unwrap', geocode=True)[0])
print (sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm, geocode=True)[0])
# Incidence angle
print (sbas.incidence_angle())
# Displacements
print (sbas.open_grids(pairs, 'unwrap', func=[sbas.vertical_displacement_mm, sbas.nearest_grid], mask=sbas.intf_ra2ll(landmask_ra), geocode=True)[0])
print (sbas.open_grids(pairs, 'unwrap', func=[sbas.eastwest_displacement_mm, sbas.nearest_grid], mask=sbas.intf_ra2ll(landmask_ra), geocode=True)[0])
