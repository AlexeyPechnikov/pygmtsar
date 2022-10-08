#!/usr/bin/env python3
# download dataset from GMTSAR official site:
#wget -c http://topex.ucsd.edu/gmtsar/tar/S1A_Stack_CPGF_T173.tar.gz
#tar xvzf S1A_Stack_CPGF_T173.tar.gz -C .

from pygmtsar import SBAS
import numpy as np
import xarray as xr
import os

#Define Parameters
BACKUPDIR    = 'backup'
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
#decimator = lambda dataarray: dataarray.coarsen({'y': 4, 'x': 4}, boundary='trim').mean()
# SRTM3 DEM resolution 3 sec, 90 m but use 60 m instead to follow NetCDF chunks size 512x512
decimator = sbas.pixel_decimator(resolution_meters=60, debug=True)
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

# test results
# https://github.com/mobigroup/gmtsar/issues/2
xr.testing.assert_allclose(
    sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm),
    sbas.los_displacement_mm(sbas.open_grids(pairs, 'unwrap'))
)
xr.testing.assert_allclose(
    sbas.open_grids(pairs, 'unwrap', geocode=True, func=sbas.los_displacement_mm),
    sbas.intf_ra2ll(sbas.los_displacement_mm(sbas.open_grids(pairs, 'unwrap')))
)
xr.testing.assert_allclose(
    sbas.open_grids(None, 'vel', geocode=True, add_subswath=False),
    sbas.intf_ra2ll(sbas.open_grids(None, 'vel', add_subswath=False))
)

# dump and restore current state
sbas.dump()
sbas = SBAS.restore(WORKDIR)
assert os.path.exists(os.path.join(WORKDIR,   'SBAS.pickle')), 'State file not found'

# backup (copy)
sbas.backup('backup', copy=True)
assert os.path.exists(os.path.join(BACKUPDIR, 'SBAS.pickle')), 'Backup state file not found'
scenes = sbas.df.datapath.apply(lambda fn: os.path.splitext(os.path.split(fn)[-1])[0]).values
for scene in scenes:
    for ext in ['.tiff', '.xml']:
        assert os.path.exists(os.path.join(BACKUPDIR,   scene+ext)), 'Backup file not found'
        assert os.path.exists(os.path.join(DATADIR,     scene+ext)), 'Original file not found'
        assert not os.path.exists(os.path.join(WORKDIR, scene+ext)), 'Excessive file found'
orbits = sbas.df.orbitpath.apply(lambda fn: os.path.split(fn)[-1]).values
for orbit in orbits:
    assert os.path.exists(os.path.join(BACKUPDIR,   orbit)), 'Backup file not found'
    assert os.path.exists(os.path.join(DATADIR,     orbit)), 'Original file not found'
    assert not os.path.exists(os.path.join(WORKDIR, orbit)), 'Excessive file found'

# backup (move only local files)
sbas.backup('backup', copy=False)
for scene in scenes:
    for ext in ['.tiff', '.xml']:
        assert os.path.exists(os.path.join(BACKUPDIR,   scene+ext)), 'Backup file not found'
        assert os.path.exists(os.path.join(DATADIR,     scene+ext)), 'Original file not found'
        assert not os.path.exists(os.path.join(WORKDIR, scene+ext)), 'Excessive file found'
for orbit in orbits:
    assert os.path.exists(os.path.join(BACKUPDIR,   orbit)), 'Backup file not found'
    assert os.path.exists(os.path.join(DATADIR,     orbit)), 'Original file not found'
    assert not os.path.exists(os.path.join(WORKDIR, orbit)), 'Excessive file found'
