# -*- coding: utf-8 -*-
"""YamchiDam.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1diVdEahWeJdzvBG7tUQJrj-13uzs6feS

## PyGMTSAR SBAS and PSI Analyses: Yamchi dam, Iran

The PyGMTSAR InSAR library, Geomed3D Geophysical Inversion Library, N-Cube 3D/4D GIS Data Visualization, among others, are my open-source projects developed in my free time. I hold a Master's degree in STEM, specializing in radio physics. In 2004, I received the first prize in the All-Russian Physics Competition for significant results in forward and inverse modeling for nonlinear optics and holography. These skills are also applicable to modeling Gravity, Magnetic, and Thermal fields, as well as satellite interferometry processing. With 20 years of experience as a data scientist and software developer, I have contributed to scientific and industrial development, working on government contracts, university projects, and with companies like LG Corp and Google Inc.

You can support my work on [Patreon](https://www.patreon.com/pechnikov), where I share updates on my projects, publications, use cases, examples, and other useful information. For research and development services and support, please visit my profile on the freelance platform [Upwork](https://www.upwork.com).

### Resources
- Google Colab Pro notebooks and articles on [Patreon](https://www.patreon.com/pechnikov),
- Google Colab notebooks on [GitHub](https://github.com),
- Docker Images on [DockerHub](https://hub.docker.com),
- Geological Models on [YouTube](https://www.youtube.com),
- VR/AR Geological Models on [GitHub](https://github.com),
- Live updates and announcements on [LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/).

© Alexey Pechnikov, 2024

$\large\color{blue}{\text{Hint: Use menu Cell} \to \text{Run All or Runtime} \to \text{Complete All or Runtime} \to \text{Run All}}$
$\large\color{blue}{\text{(depending of your localization settings) to execute the entire notebook}}$

## Load Modules to Check Environment
"""

import platform, sys, os

"""## Google Colab Installation

### Install GMTSAR
https://github.com/gmtsar/gmtsar
"""

if 'google.colab' in sys.modules:
    count = !ls /usr/local | grep GMTSAR | wc -l
    if count == ['0']:
        !export DEBIAN_FRONTEND=noninteractive
        !apt-get update > /dev/null
        !apt install -y csh autoconf gfortran \
            libtiff5-dev libhdf5-dev liblapack-dev libgmt-dev gmt-dcw gmt-gshhg gmt  > /dev/null
        # GMTSAR codes are not so good to be compiled by modern GCC
        !apt install gcc-9 > /dev/null
        !update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 10
        !update-alternatives --config gcc
        !gcc --version | head -n 1
        !rm -fr /usr/local/GMTSAR
        !git config --global advice.detachedHead false
        !cd /usr/local && git clone -q --branch master https://github.com/gmtsar/gmtsar GMTSAR
        # revert recent broken commit
        !cd /usr/local/GMTSAR && git checkout e98ebc0f4164939a4780b1534bac186924d7c998 > /dev/null
        !cd /usr/local/GMTSAR && autoconf > /dev/null
        !cd /usr/local/GMTSAR && ./configure --with-orbits-dir=/tmp > /dev/null
        !cd /usr/local/GMTSAR && make 1>/dev/null 2>/dev/null
        !cd /usr/local/GMTSAR && make install >/dev/null
        # fix for missed script, use bash instead of csh interpretator
        # note: csh messes stdout and stderr in Docker environment, it's resolved in PyGMTSAR code
        !echo '#!/bin/sh' > /usr/local/GMTSAR/bin/gmtsar_sharedir.csh
        !echo echo /usr/local/GMTSAR/share/gmtsar >> /usr/local/GMTSAR/bin/gmtsar_sharedir.csh
        !chmod a+x /usr/local/GMTSAR/bin/gmtsar_sharedir.csh
        !/usr/local/GMTSAR/bin/gmtsar_sharedir.csh
        # test one GMTSAR binary
        !/usr/local/GMTSAR/bin/make_s1a_tops 2>&1 | head -n 2

import sys
if 'google.colab' in sys.modules:
    !apt install -y xvfb > /dev/null
    !{sys.executable} -m pip install pyvista xvfbwrapper > /dev/null
    import xvfbwrapper
    display = xvfbwrapper.Xvfb(width=800, height=600)
    display.start()

"""### Define ENV Variables for Jupyter Instance"""

# Commented out IPython magic to ensure Python compatibility.
# use default GMTSAR installation path
PATH = os.environ['PATH']
if PATH.find('GMTSAR') == -1:
    PATH = os.environ['PATH'] + ':/usr/local/GMTSAR/bin/'
#     %env PATH {PATH}

"""### Install Python Modules

Maybe you need to restart your notebook, follow the instructions printing below.

The installation takes a long time on fresh Debian 10 and a short time on Google Colab
"""

!{sys.executable} --version

if 'google.colab' in sys.modules:
    #!{sys.executable} -m pip install -Uq git+https://github.com/mobigroup/gmtsar.git@pygmtsar2#subdirectory=pygmtsar
    !{sys.executable} -m pip install -q pygmtsar
from pygmtsar import __version__
__version__

"""## Load and Setup Python Modules"""

import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import json
from dask.distributed import Client
import dask
import warnings
warnings.filterwarnings('ignore')

# Commented out IPython magic to ensure Python compatibility.
# plotting modules
import pyvista as pv
# magic trick for white background
pv.set_plot_theme("document")
import panel
panel.extension('vtk')
#import seaborn as sns
#import adjustText
from contextlib import contextmanager
import matplotlib.pyplot as plt
@contextmanager
def mpl_settings(settings):
    original_settings = {k: plt.rcParams[k] for k in settings}
    plt.rcParams.update(settings)
    yield
    plt.rcParams.update(original_settings)
plt.rcParams['figure.figsize'] = [12, 4]
plt.rcParams['figure.dpi'] = 150
plt.rcParams['figure.titlesize'] = 24
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
# %matplotlib inline

# define Pandas display settings
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 100)

from pygmtsar import S1, Stack, tqdm_dask, NCubeVTK, ASF, XYZTiles
if os.path.exists('/.dockerenv') and not 'google.colab' in sys.modules:
    # use different NetCDF backend in Docker containers
    from pygmtsar import datagrid
    datagrid.netcdf_engine = 'netcdf4'

"""## Define Sentinel-1 SLC Scenes and Processing Parameters

When you need more scenes and SBAS analysis  see examples on PyGMTSAR GitHub page https://github.com/mobigroup/gmtsar

### Ascending Orbit Configuration

https://search.asf.alaska.edu/#/?polygon=POINT(48.0824%2038.0694)&start=2022-12-31T17:00:00Z&resultsLoaded=true&productTypes=SLC&granule=S1A_IW_SLC__1SDV_20230110T144512_20230110T144539_046723_0599DF_4A6F-SLC&zoom=6.528&center=47.103,33.893&end=2023-12-31T16:59:59Z&path=101-101
"""

SCENES = """
S1A_IW_SLC__1SDV_20231224T144517_20231224T144545_051798_0641BE_9896
S1A_IW_SLC__1SDV_20231212T144518_20231212T144546_051623_063BB8_5761
S1A_IW_SLC__1SDV_20231130T144518_20231130T144546_051448_06358B_D5A3
S1A_IW_SLC__1SDV_20231118T144519_20231118T144547_051273_062F91_62BA
S1A_IW_SLC__1SDV_20231106T144519_20231106T144547_051098_062982_8077
S1A_IW_SLC__1SDV_20231025T144520_20231025T144548_050923_062384_3A6A
S1A_IW_SLC__1SDV_20231013T144520_20231013T144548_050748_061D8C_31FD
S1A_IW_SLC__1SDV_20231001T144520_20231001T144548_050573_06178B_869A
S1A_IW_SLC__1SDV_20230919T144520_20230919T144548_050398_061190_5888
S1A_IW_SLC__1SDV_20230907T144519_20230907T144547_050223_060B90_A1C0
S1A_IW_SLC__1SDV_20230826T144519_20230826T144547_050048_0605A0_46AB
S1A_IW_SLC__1SDV_20230814T144518_20230814T144546_049873_05FFA4_616E
S1A_IW_SLC__1SDV_20230802T144517_20230802T144545_049698_05F9E7_960C
S1A_IW_SLC__1SDV_20230721T144517_20230721T144544_049523_05F47C_BC96
S1A_IW_SLC__1SDV_20230709T144516_20230709T144544_049348_05EF28_2658
S1A_IW_SLC__1SDV_20230627T144515_20230627T144542_049173_05E9C4_B7F4
S1A_IW_SLC__1SDV_20230615T144514_20230615T144542_048998_05E46C_CD01
S1A_IW_SLC__1SDV_20230603T144514_20230603T144542_048823_05DF18_89EE
S1A_IW_SLC__1SDV_20230522T144513_20230522T144541_048648_05D9E5_E9BA
S1A_IW_SLC__1SDV_20230510T144512_20230510T144540_048473_05D4B8_4527
S1A_IW_SLC__1SDV_20230428T144512_20230428T144540_048298_05CEF3_D42C
S1A_IW_SLC__1SDV_20230416T144511_20230416T144539_048123_05C918_63D0
S1A_IW_SLC__1SDV_20230404T144511_20230404T144539_047948_05C32C_2731
S1A_IW_SLC__1SDV_20230323T144511_20230323T144538_047773_05BD44_50B0
S1A_IW_SLC__1SDV_20230311T144511_20230311T144538_047598_05B75F_FD5F
S1A_IW_SLC__1SDV_20230227T144511_20230227T144539_047423_05B175_5102
S1A_IW_SLC__1SDV_20230215T144510_20230215T144538_047248_05AB78_36F6
S1A_IW_SLC__1SDV_20230203T144511_20230203T144539_047073_05A59C_2861
S1A_IW_SLC__1SDV_20230122T144511_20230122T144539_046898_059FC7_FBC8
"""
#S1A_IW_SLC__1SDV_20230110T144512_20230110T144539_046723_0599DF_4A6F
SCENES = list(filter(None, SCENES.split('\n')))
print (f'Scenes defined: {len(SCENES)}')

ORBIT     = 'A'
SUBSWATH  = 2
REFERENCE = '2023-07-21'

WORKDIR = 'raw_yamchi_desc'  if ORBIT == 'D' else 'raw_yamchi_asc'
DATADIR = 'data_yamchi_desc' if ORBIT == 'D' else 'data_yamchi_asc'

geojson = '''
{
  "type": "Feature",
  "geometry": {
    "type": "Point",
    "coordinates": [48.08244, 38.069423]
  },
  "properties": {"name": "Yamchi Dam, Iran"}
}
'''
POI = gpd.GeoDataFrame.from_features([json.loads(geojson)])
POI

BUFFER = 0.05
AOI = POI.buffer(BUFFER)
AOI

"""## Download and Unpack Datasets

## Enter Your ASF (Earthdata) User and Password

If the data directory is empty or doesn't exist, you'll need to download Sentinel-1 scenes from the Alaska Satellite Facility (ASF) datastore. Use your Earthdata Login credentials. If you don't have an Earthdata Login, you can create one at https://urs.earthdata.nasa.gov//users/new

You can also use pre-existing SLC scenes stored on your Google Drive, or you can copy them using a direct public link from iCloud Drive.

The credentials below are available at the time the notebook is validated. Special symbols, like underscores, are required in your ASF password.
"""

# Set these variables to None and you will be prompted to enter your username and password below.
username = 'GoogleColab2023'
password = 'GoogleColab_2023'

# download required polarization and subswaths only
asf = ASF(username, password)
asf.download(DATADIR, SCENES, SUBSWATH)

"""## Run Local Dask Cluster

Launch Dask cluster for local and distributed multicore computing. That's possible to process terabyte scale Sentinel-1 SLC datasets on Apple Air 16 GB RAM.
"""

# simple Dask initialization
if 'client' in globals():
    client.close()
client = Client()
client

"""## Init

Search recursively for measurement (.tiff) and annotation (.xml) and orbit (.EOF) files in the DATA directory. It can be directory with full unzipped scenes (.SAFE) subdirectories or just a directory with the list of pairs of required .tiff and .xml files (maybe pre-filtered for orbit, polarization and subswath to save disk space). If orbit files and DEM are missed these will be downloaded automatically below.

### Select Original Secenes and Orbits

Use filters to find required subswath, polarization and orbit in original scenes .SAFE directories in the data directory.
"""

scenes = S1.scan_slc(DATADIR, subswath=SUBSWATH)

sbas = Stack(WORKDIR, drop_if_exists=True).set_scenes(scenes).set_reference(REFERENCE)
sbas.to_dataframe()

sbas.plot_scenes(AOI=AOI)

"""## Reframe Scenes (Optional)

Stitch sequential scenes and crop the subswath to a smaller area for faster processing when the full area is not needed.
"""

sbas.compute_reframe(AOI)

sbas.plot_scenes(AOI=AOI)

"""### Download SRTM DEM

The function below downloads SRTM1 or SRTM3 DEM and converts heights to ellipsoidal model using EGM96 grid.
Besides, for faster processing we can use pre-defined DEM file as explained above. Select product=SRTM1 for 30m resolution and product=SRTM3 for 90m resolution SRTM DEM.
"""

# define the area of interest (AOI) to speedup the processing
sbas.download_dem(AOI)

sbas.plot_scenes(AOI=AOI)

"""## Align Images"""

if os.path.exists('/.dockerenv') and not 'google.colab' in sys.modules:
    # use special joblib backend in Docker containers
    sbas.compute_align(joblib_aligning_backend='threading')
else:
    sbas.compute_align()

"""## Backup"""

# bursts-cropped GeoTIFF files moved into backup directory
sbas.backup('backup')

# free disk space removing the cropped Sentinel-1 GeoTIFFs
# alternatively, drop the data directory and use the backup
!rm -fr backup

"""## Geocoding Transform"""

# use the original Sentinel-1 resolution (1 pixel spacing)
sbas.compute_geocode(1)

sbas.plot_topo(quantile=[0.01, 0.99])

"""## SBAS Baseline"""

baseline_pairs = sbas.sbas_pairs(days=200, meters=150)
# optionally, drop dates having less then 2 pairs
baseline_pairs = sbas.sbas_pairs_limit(baseline_pairs, limit=2, iterations=2)
baseline_pairs

with mpl_settings({'figure.dpi': 300}):
    sbas.plot_baseline(baseline_pairs)

"""## Persistent Scatterers Function (PSF)"""

# use the only selected dates for the pixels stability analysis
sbas.compute_ps()

sbas.plot_psfunction(quantile=[0.01, 0.90])

# export PSF grid to view in. QGIS, ArcGIS, etc.
# use .compute() to prevent lazy dataset saving issues
!rm -f 'psfunction.nc'
sbas.ra2ll(sbas.psfunction()).compute().to_netcdf('psfunction.nc', engine=sbas.netcdf_engine)

"""## Interferogram

The code below is detailed for education reasons and can be more compact excluding optional arguments. See other PyGMTSAR examples for shorter version.

### Multi-looked Resolution for SBAS
"""

sbas.compute_interferogram_multilook(baseline_pairs, 'intf_mlook', wavelength=100, weight=sbas.psfunction())

# optionally, materialize to disk and open
ds_sbas = sbas.open_stack('intf_mlook')
intf_sbas = ds_sbas.phase
corr_sbas = ds_sbas.correlation
corr_sbas

sbas.plot_interferograms(intf_sbas[:8], caption='SBAS Phase, [rad]')

sbas.plot_correlations(corr_sbas[:8], caption='SBAS Correlation')

"""## SBAS Analysis

### 2D Unwrapping
"""

corr_sbas_stack = corr_sbas.mean('pair')

corr_sbas_stack = sbas.sync_cube(corr_sbas_stack, 'corr_sbas_stack')

#CORRLIMIT = np.nanmedian(corr_sbas_stack)
CORRLIMIT = 0.4
sbas.plot_correlation_stack(corr_sbas_stack, CORRLIMIT, caption='SBAS Correlation Stack')

unwrap_sbas = sbas.unwrap_snaphu(intf_sbas, corr_sbas).where(corr_sbas_stack>=CORRLIMIT)
unwrap_sbas

# optionally, materialize to disk and open
unwrap_sbas = sbas.sync_cube(unwrap_sbas, 'unwrap_sbas')

sbas.plot_phases(unwrap_sbas.phase[:8], caption='SBAS Phase, [rad]')

"""### Trend Correction"""

decimator = sbas.decimator(resolution=15, grid=(1,1))
topo = decimator(sbas.get_topo())
inc = decimator(sbas.incidence_angle())
yy, xx = xr.broadcast(topo.y, topo.x)
trend_sbas = sbas.regression(unwrap_sbas.phase,
        [topo,    topo*yy,    topo*xx,    topo*yy*xx,
         topo**2, topo**2*yy, topo**2*xx, topo**2*yy*xx,
         topo**3, topo**3*yy, topo**3*xx, topo**3*yy*xx,
         inc,     inc**yy,    inc*xx,     inc*yy*xx,
         yy, xx,
         yy**2, xx**2, yy*xx,
         yy**3, xx**3, yy**2*xx, xx**2*yy], corr_sbas)

# optionally, materialize to disk and open
trend_sbas = sbas.sync_cube(trend_sbas, 'trend_sbas')

sbas.plot_phases(trend_sbas[:8], caption='SBAS Trend Phase, [rad]', quantile=[0.01, 0.99])

sbas.plot_phases((unwrap_sbas.phase - trend_sbas)[:8], caption='SBAS Phase - Trend, [rad]', vmin=-np.pi, vmax=np.pi)

"""### Quality Check"""

baseline_pairs['stddev'] = (unwrap_sbas.phase - trend_sbas).std(['y', 'x'])
print (len(baseline_pairs))
baseline_pairs

pairs_best = sbas.sbas_pairs_covering_deviation(baseline_pairs, 5)
print (len(pairs_best))
pairs_best

with mpl_settings({'figure.dpi': 300}):
    sbas.plot_baseline(pairs_best)

sbas.plot_baseline_deviation(baseline_pairs, pairs_best)

sbas.plot_baseline_duration(baseline_pairs, column='stddev', ascending=True)

sbas.plot_baseline_duration(pairs_best, column='stddev', ascending=True)

"""### Exclude Noisy Interferogram"""

trend_sbas_best  = trend_sbas.sel(pair=pairs_best.pair.values)
unwrap_sbas_best = unwrap_sbas.sel(pair=pairs_best.pair.values)
intf_sbas_best   = intf_sbas.sel(pair=pairs_best.pair.values)
corr_sbas_best   = corr_sbas.sel(pair=pairs_best.pair.values)

sbas.plot_phases((unwrap_sbas_best.phase - trend_sbas_best)[:8], caption='SBAS Phase - Trend, [rad]', vmin=-np.pi, vmax=np.pi)

"""### Coherence-Weighted Least-Squares Solution for LOS Displacement, mm"""

# calculate phase displacement in radians and convert to LOS displacement in millimeter
disp_sbas = sbas.los_displacement_mm(sbas.lstsq(
    unwrap_sbas_best.phase - trend_sbas_best,
    corr_sbas_best
))

# optionally, materialize to disk and open
disp_sbas = sbas.sync_cube(disp_sbas, 'disp_sbas')

sbas.plot_displacements(disp_sbas[::3], caption='SBAS Cumulative LOS Displacement, [mm]', quantile=[0.01, 0.99])

"""### STL model for LOS Displacement, mm"""

stl_sbas = sbas.stl(disp_sbas)
stl_sbas

stl_sbas = sbas.sync_cube(stl_sbas, 'stl_sbas')

years = ((stl_sbas.date.max() - stl_sbas.date.min()).dt.days/365.25).item()
print ('years', np.round(years, 3))
velocity_sbas = stl_sbas.trend.mean('date')/years
velocity_sbas

sbas.plot_velocity(sbas.as_geo(sbas.ra2ll(velocity_sbas)).rio.clip(AOI.geometry.envelope),
                   caption='SBAS LOS Velocity STL Decompose, 2021',
                   quantile=[0.01, 0.99], aspect='equal', POI=POI, marker='x')

fig = plt.figure(figsize=(12,4), dpi=300)

zmin, zmax = np.nanquantile(velocity_sbas, [0.01, 0.99])

ax = fig.add_subplot(1, 2, 1)
velocity_sbas.plot.imshow(cmap='turbo', vmin=zmin, vmax=zmax, ax=ax)
sbas.geocode(AOI.boundary).plot(ax=ax)
sbas.geocode(POI).plot(ax=ax, marker='x', c='r', markersize=100, label='POI')
ax.set_aspect('auto')
ax.legend(fontsize=14)
ax.set_title('Velocity, mm/year', fontsize=16)

ax = fig.add_subplot(1, 2, 2)
sbas.as_geo(sbas.ra2ll(velocity_sbas)).rio.clip(AOI.geometry)\
    .plot.imshow(cmap='turbo', vmin=zmin, vmax=zmax, ax=ax)
AOI.boundary.plot(ax=ax)
POI.plot(ax=ax, marker='x', c='r', markersize=100, label='POI')
ax.legend(loc='upper left', fontsize=14)
ax.set_title('Velocity, mm/year', fontsize=16)

plt.suptitle('SBAS LOS Velocity STL Decompose, 2021', fontsize=18)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 4), dpi=300)

x, y = [(geom.x.item(), geom.y.item()) for geom in sbas.geocode(POI).geometry][0]
disp_pixel = disp_sbas.sel(y=y, x=x, method='nearest')
stl_pixel = stl_sbas.sel(y=y, x=x, method='nearest')
plt.plot(disp_pixel.date, disp_pixel, c='r', lw=2, label='Displacement POI')
plt.plot(stl_pixel.date, stl_pixel.trend, c='r', ls='--', lw=2, label='Trend POI')
plt.plot(stl_pixel.date, stl_pixel.seasonal, c='r', lw=1, label='Seasonal POI')

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=14)
plt.title('SBAS LOS Displacement STL Decompose, 2021', fontsize=18)
plt.ylabel('Displacement, mm', fontsize=16)
plt.show()

"""## PS Analysis

Use the trend detected on possibly lower resolution unwrapped phases for higher resolution analysis.
"""

sbas.compute_interferogram_singlelook(pairs_best, 'intf_slook', wavelength=100,
                                      weight=sbas.psfunction(), phase=trend_sbas_best)

# optionally, materialize to disk and open
ds_ps = sbas.open_stack('intf_slook')
intf_ps = ds_ps.phase
corr_ps = ds_ps.correlation

"""### Exclude Noisy Interferogram

intf_ps = intf_ps.sel(pair=pairs_best.pair.values)
corr_ps = corr_ps.sel(pair=pairs_best.pair.values)

### 1D Unwrapping for LOS Displacement, mm
"""

corr_ps_stack = corr_ps.mean('pair')

corr_ps_stack = sbas.sync_cube(corr_ps_stack, 'corr_ps_stack')

sbas.plot_correlation_stack(corr_ps_stack, 0, caption='PS Stack Correlation')

disp_ps_pairs = sbas.los_displacement_mm(sbas.unwrap1d(intf_ps))
disp_ps_pairs

disp_ps_pairs = sbas.sync_cube(disp_ps_pairs, 'disp_ps_pairs')

"""### Coherence-Weighted Least-Squares Solution for LOS Displacement, mm"""

disp_ps = sbas.lstsq(disp_ps_pairs, corr_ps)
disp_ps

disp_ps = sbas.sync_cube(disp_ps, 'disp_ps')

zmin, zmax = np.nanquantile(disp_ps, [0.01, 0.99])
sbas.plot_displacements(disp_ps[::3], caption='PS Cumulative LOS Displacement, [mm]',vmin=zmin, vmax=zmax)

"""### STL model for LOS Displacement, mm"""

stl_ps = sbas.stl(disp_ps)
stl_ps

stl_ps = sbas.sync_cube(stl_ps, 'stl_ps')

years = ((stl_ps.date.max() - stl_ps.date.min()).dt.days/365.25).item()
print ('years', np.round(years, 3))
velocity_ps = stl_ps.trend.mean('date')/years
velocity_ps

sbas.plot_velocity(sbas.as_geo(sbas.ra2ll(velocity_ps)).rio.clip(AOI.geometry.envelope),
                   caption='PS LOS Velocity STL Decompose, 2021',
                   quantile=[0.01, 0.99], aspect='equal', POI=POI, marker='x')

x, y = [(geom.x.item(), geom.y.item()) for geom in sbas.geocode(POI).geometry][0]
sbas.plot_baseline_displacement(disp_ps_pairs.sel(y=y, x=x, method='nearest')/sbas.los_displacement_mm(1),
                                corr_ps.sel(y=y, x=x, method='nearest'),
                               caption='POI', stl=True)

rmse = sbas.rmse(disp_ps_pairs, disp_ps, corr_ps)
rmse

sbas.plot_rmse(sbas.as_geo(sbas.ra2ll(rmse)).rio.clip(AOI.geometry.envelope),
                   quantile=[0.01, 0.99], aspect='equal', POI=POI, marker='x')

rmselimit = 0.15
sbas.plot_velocity(sbas.as_geo(sbas.ra2ll(velocity_ps.where(rmse<rmselimit))).rio.clip(AOI.geometry.envelope),
                   caption=f'PS LOS Velocity STL Decompose RMSE<{rmselimit}',
                   quantile=[0.01, 0.99], aspect='equal', POI=POI, marker='x')

"""## SBAS vs PS Comparision"""

# crop AOI
points_sbas = sbas.as_geo(sbas.ra2ll(velocity_sbas)).rio.clip(AOI.geometry)
points_ps = sbas.as_geo(sbas.ra2ll(velocity_ps.where(rmse<rmselimit))).rio.clip(AOI.geometry)
points_ps = points_ps.interp_like(points_sbas, method='nearest').values.ravel()
points_sbas = points_sbas.values.ravel()
nanmask = np.isnan(points_sbas) | np.isnan(points_ps)
points_sbas = points_sbas[~nanmask]
points_ps = points_ps[~nanmask]

plt.figure(figsize=(12, 4), dpi=300)
plt.scatter(points_sbas, points_ps, c='silver', alpha=1,   s=1)
plt.scatter(points_sbas, points_ps, c='b',      alpha=0.1, s=1)
plt.scatter(points_sbas, points_ps, c='g',      alpha=0.1, s=0.5)
plt.scatter(points_sbas, points_ps, c='y',      alpha=0.1, s=0.25)

# adding a 1:1 line
max_value = max(velocity_sbas.max(), velocity_ps.max())
min_value = min(velocity_sbas.min(), velocity_ps.min())
plt.plot([min_value, max_value], [min_value, max_value], 'k--')

plt.xlabel('Velocity SBAS, mm/year', fontsize=16)
plt.ylabel('Velocity PS, mm/years', fontsize=16)
plt.title('Cross-Comparison between SBAS and PS Velocity', fontsize=18)
plt.grid(True)
plt.show()

"""## 3D Interactive Map"""

dem = sbas.get_dem()

velocity_sbas_ll = sbas.ra2ll(velocity_sbas)
velocity_ps_ll = sbas.ra2ll(velocity_ps.where(rmse<rmselimit))

velocity_sbas_ll = sbas.as_geo(velocity_sbas_ll).rio.clip(AOI.geometry.envelope)
velocity_ps_ll = sbas.as_geo(velocity_ps_ll).rio.clip(AOI.geometry.envelope)

gmap_tiles = XYZTiles().download(velocity_sbas_ll, 14)

for name, velocity_ll in {'sbas': velocity_sbas_ll, 'ps': velocity_ps_ll}.items():
    #sbas.as_geo(velocity_ll).rio.clip(AOI.geometry.buffer(-BUFFER))
    gmap = gmap_tiles.interp_like(velocity_ll, method='cubic').round().astype(np.uint8)

    #zmin, zmax = np.nanquantile(velocity_ll, [0.01, 0.99])
    #velocity_ll_norm = np.clip((velocity_ll - zmin) / (zmax - zmin), 0, 1)
    #cmap = plt.cm.get_cmap('turbo')
    #velocity_ll_colors = xr.DataArray(255*cmap(velocity_ll_norm)[..., :3].transpose(2,0,1), coords=gmap.coords)\
    #    .round().astype(np.uint8)\
    #    .rename('colors')
    #mask_invert = sbas.as_geo(velocity_ll_colors).rio.clip(AOI.geometry, invert=True)
    #map_colors = xr.where(mask_invert.sum('band') == 0, velocity_ll_colors, gmap)
    ds = xr.merge([dem.interp_like(velocity_ll, method='cubic').rename('z'), gmap, velocity_ll])
    # decimate large grid
    ds = ds.sel(lon=ds.lon[::5])
    # convert to VTK structure
    vtk_grid = pv.StructuredGrid(NCubeVTK.ImageOnTopography(ds.rename({'lat': 'y', 'lon': 'x'})))
    vtk_grid.save(f'velocity_{name}.vtk')
vtk_grid

plotter = pv.Plotter(shape=(1, 2), notebook=True)
axes = pv.Axes(show_actor=True, actor_scale=2.0, line_width=5)

plotter.subplot(0, 0)
vtk_grid = pv.read('velocity_sbas.vtk')
mesh = vtk_grid.scale([1, 1, 0.00002]).rotate_z(135, point=axes.origin)
plotter.add_mesh(mesh.scale([1, 1, 0.999]), scalars='colors', rgb=True, ambient=0.2)
plotter.add_mesh(mesh, scalars='trend', ambient=0.2, cmap='turbo', clim=(-20,20), nan_opacity=0.1, nan_color='black')
plotter.show_axes()
plotter.add_title('SBAS LOS Velocity', font_size=32)

plotter.subplot(0, 1)
vtk_grid = pv.read('velocity_ps.vtk')
mesh = vtk_grid.scale([1, 1, 0.00002]).rotate_z(135, point=axes.origin)
plotter.add_mesh(mesh.scale([1, 1, 0.999]), scalars='colors', rgb=True, ambient=0.2)
plotter.add_mesh(mesh, scalars='trend', ambient=0.2, cmap='turbo', clim=(-20,20), nan_opacity=0.1, nan_color='black')
plotter.show_axes()
plotter.add_title('PS LOS Velocity', font_size=32)

plotter.show_axes()
plotter._on_first_render_request()
panel.panel(
    plotter.render_window, orientation_widget=plotter.renderer.axes_enabled,
    enable_keybindings=False, sizing_mode='stretch_width', min_height=600
)

"""## Export VTK file from Google Colab"""

if 'google.colab' in sys.modules:
    from google.colab import files
    files.download('velocity_sbas.vtk')
    files.download('velocity_ps.vtk')

"""## Conclusion

For now you have the full control on interferometry processing and unwrapping and able to run it everywhere: on free of charge Google Colab instances, on local MacOS and Linux computers and on Amazon EC2 and Google Cloud VM and AI Notebook instances.
"""