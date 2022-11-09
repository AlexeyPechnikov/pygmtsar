#!/usr/bin/env python3
# S1A_IW_SLC__1SSV_20160408T091355_20160408T091430_010728_01001F_83EB.zip
# S1A_IW_SLC__1SSV_20160420T091355_20160420T091423_010903_010569_F9CE.zip

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
    #decimator = lambda dataarray: dataarray.coarsen({'y': 4, 'x': 4}, boundary='trim').mean()
    # use the DEM default resolution 60 m
    decimator = sbas.pixel_decimator(resolution_meters=60, debug=True)
    sbas.intf_parallel(pairs, func=decimator)
    # geocode matrices
    sbas.geocode_parallel(pairs)
    # output
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
    sbas.sat_look_parallel()
    print (sbas.incidence_angle())
    # Displacements
    print (sbas.open_grids(pairs, 'unwrap', func=[sbas.vertical_displacement_mm, sbas.nearest_grid], mask=landmask_ra, geocode=True)[0])
    print (sbas.open_grids(pairs, 'unwrap', func=[sbas.eastwest_displacement_mm, sbas.nearest_grid], mask=landmask_ra, geocode=True)[0])

    phasefilt = sbas.open_grids(pairs, 'phasefilt', geocode=True)[0]
    plt.figure(figsize=(12,4), dpi=300)
    phasefilt.plot.imshow(vmin=-np.pi, vmax=np.pi, cmap='gist_rainbow_r')
    plt.title('Phase, [rad]', fontsize=18)
    plt.savefig('Phase, [rad].jpg', dpi=150, pil_kwargs={"quality": 95})

    vert_disp_mm = sbas.open_grids(pairs, 'unwrap', func=[sbas.vertical_displacement_mm, sbas.nearest_grid], mask=sbas.intf_ra2ll(landmask_ra), geocode=True)[0]
    plt.figure(figsize=(12,4), dpi=300)
    zmin, zmax = np.nanquantile(vert_disp_mm, [0.01, 0.99])
    vert_disp_mm.plot.imshow(vmin=zmin, vmax=zmax, cmap='jet')
    plt.title('Landmasked Vertical Displacement, [mm]', fontsize=18)
    plt.savefig('Landmasked Vertical Displacement, [mm].jpg', dpi=150, pil_kwargs={"quality": 95})

    east_disp_mm = sbas.open_grids(pairs, 'unwrap', func=[sbas.eastwest_displacement_mm, sbas.nearest_grid], mask=sbas.intf_ra2ll(landmask_ra), geocode=True)[0]
    plt.figure(figsize=(12,4), dpi=300)
    zmin, zmax = np.nanquantile(east_disp_mm, [0.01, 0.99])
    east_disp_mm.plot.imshow(vmin=zmin, vmax=zmax, cmap='jet')
    plt.title('Landmasked East-West Displacement, [mm]', fontsize=18)
    plt.savefig('Landmasked East-West Displacement, [mm].jpg', dpi=150, pil_kwargs={"quality": 95})
