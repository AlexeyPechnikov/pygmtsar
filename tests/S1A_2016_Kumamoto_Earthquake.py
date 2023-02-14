#!/usr/bin/env python3
# S1A_IW_SLC__1SSV_20160408T091355_20160408T091430_010728_01001F_83EB.zip
# S1A_IW_SLC__1SSV_20160420T091355_20160420T091423_010903_010569_F9CE.zip

from pygmtsar import SBAS
import numpy as np
import xarray as xr
from dask.distributed import Client
import matplotlib.pyplot as plt

#Define Parameters
WORKDIR      = 'raw_kumamoto'
DATADIR      = 'data_kumamoto'
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
    # commented code below is not safe
    #print (sbas.open_grids(pairs, 'unwrap', func=[sbas.vertical_displacement_mm, sbas.nearest_grid], mask=landmask_ra, geocode=True)[0])
    #print (sbas.open_grids(pairs, 'unwrap', func=[sbas.eastwest_displacement_mm, sbas.nearest_grid], mask=landmask_ra, geocode=True)[0])
    # fill small gaps, georeference grid, calculate displacements, mask using landmask in geographic coordinates
    interp = lambda grid: sbas.nearest_grid(grid, search_radius_pixels=10)
    # landmask in geographic coordinates applied after all the user functions (when the grid is geocoded)
    landmask = sbas.get_landmask()
    print (sbas.open_grids(pairs, 'unwrap', func=[interp, sbas.intf_ra2ll, sbas.vertical_displacement_mm], mask=landmask)[0])
    print (sbas.open_grids(pairs, 'unwrap', func=[interp, sbas.intf_ra2ll, sbas.eastwest_displacement_mm], mask=landmask)[0])

    topo_ra = sbas.get_topo_ra()
    plt.figure(figsize=(12,4), dpi=300)
    topo_ra.plot.imshow(cmap='gray', vmin=0)
    plt.xlabel('Range', fontsize=16)
    plt.ylabel('Azimuth', fontsize=16)
    plt.title('Topography in Radar Coordinates, [m]', fontsize=18)
    plt.savefig('Topography in Radar Coordinates, [m].jpg', dpi=300, pil_kwargs={"quality": 95})

    phasefilt = sbas.open_grids(pairs, 'phasefilt', geocode=True)[0]
    plt.figure(figsize=(12,4), dpi=300)
    phasefilt.plot.imshow(vmin=-np.pi, vmax=np.pi, cmap='gist_rainbow_r')
    plt.title('Phase, [rad]', fontsize=18)
    plt.savefig('Phase, [rad].jpg', dpi=300, pil_kwargs={"quality": 95})

    # landmask in geographic coordinates
    landmask = sbas.get_landmask()
    plt.figure(figsize=(12,4), dpi=300)
    landmask.plot.imshow(vmin=0, vmax=1, cmap='gray')
    plt.title('Landmask', fontsize=18)
    plt.savefig('Landmask.jpg', dpi=300, pil_kwargs={"quality": 95})

    # convert landmask to radar coordinates on the fly
    landmask_ra = sbas.get_landmask(inverse_geocode=True)
    plt.figure(figsize=(12,4), dpi=300)
    landmask_ra.plot.imshow(vmin=0, vmax=1, cmap='gray')
    plt.title('Landmask in Radar Coordinates', fontsize=18)
    plt.savefig('Landmask in Radar Coordinates.jpg', dpi=300, pil_kwargs={"quality": 95})

    # landmask in geographic coordinates applied after all the user functions (when the grid is geocoded)
    landmask = sbas.get_landmask()
    phasefilt = sbas.open_grids(pairs, 'phasefilt', mask=landmask, geocode=True)[0]
    plt.figure(figsize=(12,4), dpi=300)
    phasefilt.plot.imshow(vmin=-np.pi, vmax=np.pi, cmap='gist_rainbow_r')
    plt.title('Landmasked Phase, [rad]', fontsize=18)
    plt.savefig('Landmasked Phase, [rad].jpg', dpi=300, pil_kwargs={"quality": 95})

    unwrap = sbas.open_grids(pairs, 'unwrap', geocode=True)[0]
    plt.figure(figsize=(12,4), dpi=300)
    unwrap.plot.imshow(cmap='jet')
    plt.title('Landmasked Unwrapped Phase, [rad]', fontsize=18)
    plt.savefig('Landmasked Unwrapped Phase, [rad].jpg', dpi=300, pil_kwargs={"quality": 95})

    los_disp_mm = sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm, geocode=True)[0]
    plt.figure(figsize=(12,4), dpi=300)
    zmin, zmax = np.nanquantile(los_disp_mm, [0.01, 0.99])
    los_disp_mm.plot.imshow(vmin=zmin, vmax=zmax, cmap='jet')
    plt.title('Landmasked LOS Displacement, [mm]', fontsize=18)
    plt.savefig('Landmasked LOS Displacement, [mm].jpg', dpi=300, pil_kwargs={"quality": 95})

    incidence_angle = sbas.incidence_angle()
    plt.figure(figsize=(12,4), dpi=300)
    incidence_angle.plot.imshow(cmap='gray')
    plt.title('Incidence Angle, [rad]', fontsize=18)
    plt.savefig('Incidence Angle, [rad].jpg', dpi=300, pil_kwargs={"quality": 95})

    # fill small gaps, georeference grid, calculate displacements, mask using landmask in geographic coordinates
    interp = lambda grid: sbas.nearest_grid(grid, search_radius_pixels=10)
    # landmask in geographic coordinates applied after all the user functions (when the grid is geocoded)
    landmask = sbas.get_landmask()
    vert_disp_mm = sbas.open_grids(pairs, 'unwrap', func=[interp, sbas.intf_ra2ll, sbas.vertical_displacement_mm], mask=landmask)[0]
    plt.figure(figsize=(12,4), dpi=300)
    zmin, zmax = np.nanquantile(vert_disp_mm, [0.01, 0.99])
    vert_disp_mm.plot.imshow(vmin=zmin, vmax=zmax, cmap='jet')
    plt.title('Landmasked Vertical Displacement, [mm]', fontsize=18)
    plt.savefig('Landmasked Vertical Displacement, [mm].jpg', dpi=300, pil_kwargs={"quality": 95})

    # fill small gaps, georeference grid, calculate displacements, mask using landmask in geographic coordinates
    interp = lambda grid: sbas.nearest_grid(grid, search_radius_pixels=10)
    # landmask in geographic coordinates applied after all the user functions (when the grid is geocoded)
    landmask = sbas.get_landmask()
    east_disp_mm = sbas.open_grids(pairs, 'unwrap', func=[interp, sbas.intf_ra2ll, sbas.eastwest_displacement_mm], mask=landmask)[0]
    plt.figure(figsize=(12,4), dpi=300)
    zmin, zmax = np.nanquantile(east_disp_mm, [0.01, 0.99])
    east_disp_mm.plot.imshow(vmin=zmin, vmax=zmax, cmap='jet')
    plt.title('Landmasked East-West Displacement, [mm]', fontsize=18)
    plt.savefig('Landmasked East-West Displacement, [mm].jpg', dpi=300, pil_kwargs={"quality": 95})