#!/usr/bin/env python3
# download dataset from GMTSAR official site:
#wget -c http://topex.ucsd.edu/gmtsar/tar/S1A_Stack_CPGF_T173.tar.gz
#tar xvzf S1A_Stack_CPGF_T173.tar.gz -C .

from pygmtsar import SBAS
import numpy as np

#Define Parameters
MASTER       = '2015-04-03'
WORKDIR      = 'raw'
DATADIR      = 'raw_orig'
DEMFILE      = 'topo/dem.grd'
BASEDAYS     = 100
BASEMETERS   = 150
DEFOMAX      = 0
CORRLIMIT    = 0.075

# Init SBAS
sbas = SBAS(DATADIR, DEMFILE, WORKDIR).set_master(MASTER)
# output
print (sbas.to_dataframe())
# Align a Stack of Images
sbas.stack_parallel()
# SBAS Baseline
baseline_pairs = sbas.baseline_pairs(days=BASEDAYS, meters=BASEMETERS)
# output
print (baseline_pairs)
# DEM in Radar Coordinates
sbas.topo_ra_parallel()
# output
print (sbas.get_topo_ra())
#Interferograms
pairs = baseline_pairs[['ref_date', 'rep_date']]
decimator = lambda dataarray: dataarray.coarsen({'y': 4, 'x': 4}, boundary='trim').mean()
sbas.intf_parallel(pairs, wavelength=400, func=decimator)
# output
print (pairs)
print (sbas.open_grids(pairs, 'phasefilt'))
print (sbas.open_grids(pairs, 'phasefilt', geocode=True))
print (sbas.open_grids(pairs, 'corr'))
# Unwrapping
interpolator = lambda corr, unwrap: sbas.nearest_grid(unwrap).where(corr>=0)
sbas.unwrap_parallel(pairs, threshold=CORRLIMIT, func=interpolator)
# output
print (sbas.open_grids(pairs, 'unwrap'))
print (sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm))
print (sbas.open_grids(pairs, 'unwrap', geocode=True, func=sbas.los_displacement_mm))
# GMTSAR SBAS Displacement and Velocity
sbas.sbas(baseline_pairs, smooth=1)
# output
phasefilts = sbas.open_grids(pairs, 'phasefilt')
print (sbas.open_grids(sbas.df.index, 'disp', func=sbas.nearest_grid, mask=phasefilts.min('pair'), add_subswath=False))
phasefilts_ll = sbas.open_grids(pairs, 'phasefilt', geocode=True)
print (sbas.open_grids(sbas.df.index, 'disp', geocode=True, func=sbas.nearest_grid, mask=phasefilts_ll.min('pair'), add_subswath=False))
print (sbas.open_grids(None, 'vel', add_subswath=False))
print (sbas.open_grids(None, 'vel', geocode=True, add_subswath=False))
print (sbas.open_grids(None, 'vel', func=sbas.nearest_grid, mask=phasefilts.min('pair'), add_subswath=False))
print (sbas.open_grids(None, 'vel', geocode=True, func=sbas.nearest_grid, mask=phasefilts_ll.min('pair'), add_subswath=False))
# Detrending
sbas.detrend_parallel(pairs, wavelength=12000)
# output
print (sbas.open_grids(pairs, 'detrend'))
# PyGMTSAR SBAS Displacement
sbas.sbas_parallel(pairs)
# output
valid_area = sbas.open_grids(pairs, 'corr').min('pair')
def interp(da):
    return sbas.nearest_grid(da).where(~np.isnan(valid_area))
dates = np.unique(pairs.values.flatten())
disps = sbas.open_grids(dates, 'disp', func=[sbas.los_displacement_mm, interp])
print (sbas.cropna(sbas.intf_ra2ll(disps.cumsum('date').sel(x=slice(9000,12000), y=slice(1800,2700)))))
