## About

PyGMTSAR (Python GMTSAR) is an open source project and Python package that provides easy and fast Sentinel-1 Satellite Interferometry for everyone! While it's pure Python package under the hood it uses GMTSAR binary tools which should be installed.

Initially, PyGMTSAR is forked from GMTSAR GitHub repository and lot's of changes are made to seamlessly call all the binary tools from Python API. For now, all the changes developed for PyGMTSAR project merged into GMTSAR and the both projects are binary compatible. To prevent confusing, below PyGMTSAR means the Python package only and GMTSAR means the binary core tools (while these are available in PyGMTSAR GitHub repository too). 

The goal of the project is easy and fast satellite interferometry (InSAR) interactive and batch processing in Python scripts and Jupyter Notebooks for Sentinel-1 radar scenes everywhere as on localhost as on cloud environments like to Google Cloud VM and AI Notebooks and Amazon EC2 and even on free of charge cloud environment Google Colab. 

## Why PyGMTSAR?

PyGMTSAR itself combines powerful Python instrumentary for sophisticated multidementional processing (xarray library) and lazy calculations (dask library) plus parallel computing (dask and joblib libraries) to perform fast and interactive processing on huge datasets. And the best algorithms and numerical computation approaches applied for all the processing steps. There are progressbars and preview plots for the every step and that's easy to save intermediate results and continue work later on the same or other host. And (thanks to joblib library) that's safe to interrupt the execution at any time without memory leaks (common for dask-based solutions).

Thanks to all the powerful Python libraries and the best used algorithms PyGMTSAR is really fast and that's possible to complete SBAS analysis for 5 years on 800 interferograms in just one day even on Apple Air or Apple iMac (8 cores and 16 GB RAM) using 2 TB raw Sentinel-1 scenes. See the live Google Colab notebooks to find how dramatically PyGMTSAR enhaces the results in comparision to GMTSAR output for the same case.

And PyGMTSAR is user-friendly providing functions to download the required satellite orbit files and DEM and so on. 

Although PyGMTSAR and GMTSAR are binary compatible PyGMTSAR build rich Python API for interactive and batch computations and GMTSAR provides a set of shell scripts for the batch processing only. By the way, for development purposes PyGMTSAR includes shell scripts pipeline for SBAS processing with advantages like to cloud hosts init script and Telegram messenger processing notifications and you are able to try it in case you are fun of shell scripting (that's a joke of course but who knows...).

## Documentation

Satellite interferometry is a complex field of knowledge and usually that's requires a lot of time to understood it enough to build the meaningful results. Really, there is a big gap between GMTSAR github cloning and, for an example, valid SBAS + PSI results computation. PyGMTSAR resolves the software-related complexity providing high-level functions running in Live Jupyter notebooks on Google Colab so you have immediate access to the ready to use and working interactive examples which can be easily modified to produce your own results without hassle with software installation and configuration. You'd start from the self-documented Live Google Colab examples and view PyGMTSAR functions documentation and GMTSAR installation instructions  later.

## Tutorials: Live Examples on Google Colab

Click on the links below to run the processing in your own browser without any software installation. That's like to magic and it works.

* [ASF Downloading 2017 Iran–Iraq Earthquake vs GMTSAR GAMMA SNAP Co-Seismic Interferogram](https://colab.research.google.com/drive/12LJqlZNBUmvLlRl98rRFCbKveVPg9Ami?usp=sharing). The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** and **compares the results to GMTSAR, SNAP and GAMMA Software**. Note: replace the scene names to produce an **interferogram** and **LOS displacement** for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/177748605-788889e5-9afd-44d8-bc3c-dc6efe920ea0.png" width="50%">

* [Live Example S1A_2016_Kumamoto Earthquake_Co-Seismic Interferogram vs ESA Sentinel 1 Toolbox on Alaska Satellite Facility](https://colab.research.google.com/drive/1PyYcxvuyzhh-g4NQEbKjcfTDQhREZInn?usp=sharing). This is a single subswath processing with **landmask** applied to **interferogram**, **unwapped phase**, and **LOS, east-west, vertical displacement** results.

<img src="https://user-images.githubusercontent.com/7342379/183805898-d7c1ad76-822e-428e-9259-f19cc9e7540e.jpg" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183816622-1dacce7e-6a2f-46b9-8e67-d701f55bdd30.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183649417-7fcb7f3f-8c8d-45e8-a2c9-9293498ebada.png" width="50%">

* [Live Example S1AB 2021 Crete Earthquake Co-Seismic Interferogram vs Centre of EO Research & Satellite Remote Sensing, Greece Report](https://colab.research.google.com/drive/1ZTPV4HY-UoLvDYVx0UGh_Z3B12scSh9E?usp=sharing) This is a single **cropped subswath** processing with **landmask** applied to **interferogram**, **unwapped phase**, and **LOS, east-west, vertical displacement** results.

<img src="https://user-images.githubusercontent.com/7342379/177004287-cdd4351c-0834-42ae-8e46-9da5e8b124bf.jpg" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183645260-f8529ff3-b014-499e-ba2f-ebea4937b2c2.png" width="50%">

* [GMTSAR example dataset S1A_Stack_CPGF_T173](https://colab.research.google.com/drive/1sljxm2jAMGXynq4EYam6Siz8OLcPLN0h?usp=sharing) This example illustrates **SBAS** and **PSI** analyses and **detrending** approach to remove **atmospheric noise** to produce much better results.

<img src="https://user-images.githubusercontent.com/7342379/135814732-aa0eb142-ae54-4a57-b271-c33b5174a28e.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/189961167-bf3901e5-417c-41ce-a5ca-d1c74c239a04.png" width="50%">

* [ASF Downloading 2020 Ardabil, Iran Earthquake Co-Seismic Interferogram and LOS Displacement](https://colab.research.google.com/drive/1ZBVwlkiXMhSDS96oojpWrzTyRFIxv8Rp?usp=sharing). The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **crop the area** and **merge subswaths** and **detrend** results. Note: replace the scene names to produce an interferogram for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/190452115-fd703b7f-1e9b-49e5-9d17-8d9fae563aae.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/190451656-386d6cb8-f536-447c-8274-71d4f0435408.png" width="50%">



## Reference

PyGMTSAR defines a high-level Python SBAS class for an interferometry processing. Below listed the SBAS user-friendly functions grouped by common tasks in order of the full pipeline processing steps. For simple pipelines some of the steps can be missed.

Note: not listed internal function are commonly sonamed to the related wrapped or replaced GMTSAR binaries and their documentation available interactively in docstring format for developers only.

### SBAS class functions

Firstly, import the class from PyGMTSAR Python package as

```
from pygmtsar import SBAS
```

#### Initialisation (PyGMTSAR original)

```
def __init__(self, datadir, dem_filename=None, basedir='.', landmask_filename=None,
            filter_orbit=None, filter_mission=None, filter_subswath=None, filter_polarization=None, force=True):

    Initialise SBAS object using scenes directory and some filters.

    Args:
        datadir: An directory where stores raw Sentinel-1 scenes and optionally orbits.
        dem_filename: optional WGS84 DEM NetCDF file name.
        basedir: processing directory where all the temporary and output files will be created.
        landmask_filename: optional land mask NetCDF file name.
        filter_orbit: optional filter for ascending ('A') or descending ('D') orbit.
        filter_mission: optional mission name filter ('S1A' or 'S1B' for now).
        filter_subswath: optional subswath filter (anyone from the list 1,2,3,12,23,123).
        filter_polarization: optional polarization filter (anyone from list 'VV','VH','HH','HV').
        force: optional boolean flag to drop the processing directory if exists and create an empty one.

    Returns:
        Initialised SBAS object.

    Examples:
        sbas = SBAS('data', basedir='raw')
        sbas = SBAS('data', 'data/DEM_WGS84.nc', 'raw')
```

```
def set_master(self, master):

    Define master scene for SBAS object. When this call is missed then the first scene will be master.

    Args:
        master: date string.
        
    Returns:
        SBAS object.
        
    Examples:
    		sbas = SBAS(...).set_master('2022-01-20')
```

#### Manage orbits (PyGMTSAR original)

```
def download_orbits(self):

    Download missed orbits for all the SBAS scenes.

    Args:
        None

    Returns:
        None

    Examples:
        sbas.download_orbits()
```

#### List scenes and orbits (PyGMTSAR original)

```
def to_dataframe(self):

    Return Pandas Dataframe for all SBAS scenes.

    Args:
        None
        
    Returns:
        Pandas Dataframe.
        
    Examples:
    		df = sbas.to_dataframe()
```

#### Manage DEM (PyGMTSAR original)

```
def set_dem(self, dem_filename):

    Define existing WGS84 DEM file in geographic coordinates for SBAS object.

    Args:
        dem_filename: WGS84 NetCDF DEM file in geographic coordinates.
        
    Returns:
        SBAS object.
        
    Examples:
    		sbas = sbas.set_dem('data/DEM_WGS84.nc')
    		Also, the same result is possible on SBAS initialisation:
    		sbas = SBAS(..., dem_filename='data/DEM_WGS84.nc')
```

```
def download_dem(self, backend=None, product='SRTM1', resolution_meters=60, method='cubic',
								 buffer_degrees=0.02, debug=False):

    Download missed DEM file in geographic coordinates and convert it for WGS84 ellipsoid.

    Args:
        backend: optional argument to define processing backend. When it's not defined custom
        downloading method used (based on Python "elevation" libarary) available even on old systems
        but it's slow. 'GMT' backend is much faster but requires modern GMT package installation.
        Use backend=None for old systems like to Google Colab Ubuntu 18.04 (bionic).
        product: is the DEM product from two available ones 'SRTM1' or 'SRTM3'.
        resolution_meters: the DEM grid resolution in meters. The same grid is ised for geocoded
        results output.
        method: interpolation method applied when the original DEM is lower than defined
        "resolution_meters". 'cubic' interpolation prefered because it saves well the DEM
        spatial spectrum.
        buffer_degrees: add the buffer area around the master scene approximate extent.
        debug: boolean flag to print debug information.
        
    Returns:
        None
        
    Examples:
    		Download 30m STRM1 DEM on modern systems and convert it to default 60m grid:
    		sbas.download_dem(backend='GMT')
    		
    		Download 30m STRM1 DEM on modern systems:
    		sbas.download_dem(backend='GMT', resolution_meters=30)
    		
    		Download 30m STRM1 DEM and convert it to default 60m grid on old Google Colab Ubuntu 18.04 (bionic):
    		sbas.download_dem()
    		
    		Download 90m STRM3 DEM on old Google Colab Ubuntu 18.04 (bionic):
    		sbas.download_dem('STRM3', resolution_meters=90)
```

```
def get_dem(self, subswath=None, geoloc=False, buffer_degrees=0.02):

		Return WGS84 DEM in geographic coordinates as just one or list of Xarray Daraarray.
		
    Args:
    		subswath - optional subswath number. 
        geoloc: boolen flag to crop the DEM using the master scene approximate extent.
        buffer_degrees: when geoloc=True add some buffer area around the master scene approximate extent.

    Returns:
        2D Xarray Daraarray object or list of 2D Xarray Daraarray object.

    Examples:
    		Get DEM for all the processed subswaths:
        topo_ll = sbas.get_dem()
        
        Get DEM for a single subswath IW1:
        topo_ll = sbas.get_dem(1)
```

#### Manage topography in radar coordinates (combined GMTSAR wrapper and PyGMTSAR original)

```
def topo_ra_parallel(self, n_jobs=-1, method='cubic'):

		Build WGS84 DEM in radar coordinates for interferogram processing.
		
    Args:
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
        method: interpolation method applied when the original DEM is lower than defined
        "resolution_meters". 'cubic' interpolation prefered because it saves well the DEM
        spatial spectrum.

		Returns:
        None

    Examples:
    		sbas.topo_ra_parallel()
```

```
def get_topo_ra(self, subswath=None):

		Return WGS84 DEM in radar coordinates as just one or list of Xarray Daraarray.
		
    Args:
    		subswath - optional subswath number.

    Returns:
        2D Xarray Daraarray object or list of 2D Xarray Daraarray object.

    Examples:
    		Get DEM for all the processed subswaths:
        topo_ra = sbas.get_topo_ra()
        
        Get DEM for single subswath IW1:
        topo_ra = sbas.get_topo_ra(1)
```

#### Framing (combined GMTSAR wrapper and PyGMTSAR original)

```
def set_pins(self, *args):

    Define pins for scene framing. Pins are lon/lat points to crop extra bursts minimizing processing area.

    Args:
        *args - list of pairs of lon/lat coordinate points.

    Returns:
        None

    Examples:
        Set pins for a single subswath IW1 processing:
        sbas.set_pins([47.5, 37.8, 47.5, 38.3])
        
    		Set pins for a single subswath IW2 processing:
        sbas.set_pins([48.5, 38, 48.5, 38.5])
        
        Set pins for two subswaths IW1 and IW2 processing:
        sbas.set_pins([47.5, 37.8, 47.5, 38.3], [48.5, 38, 48.5, 38.5])
        
        Set pins for all three subswaths IW1, IW2, IW3 processing:
        sbas.set_pins([47.5, 38-0.2, 47.5, 39-0.2], [48.5, 38, 48.5, 39], [49.5, 38+0.2, 49.5, 39+0.2])
```

```
def get_pins(self, subswath=None):

    Return pins list. Pins are lon/lat points to crop extra bursts minimizing processing area.

    Args:
        subswath - optional subswath number.

    Returns:
        Python list of pairs of lon/lat coordinate points.

    Examples:
        Get pins for a single subswath IW1:
        sbas.get_pins(1)
        >>> [47.5, 37.8, 47.5, 38.3]
        
    		Get pins for a single subswath IW2:
        sbas.get_pins(2)
        >>> [48.5, 38, 48.5, 38.5]
        
        Get pins for all the subswaths for a single subswath IW1 processing:
        sbas.get_pins()
        >>> [47.5, 37.8, 47.5, 38.3]
        
        Get pins for all the subswaths for a single subswath IW2 processing:
        sbas.get_pins()
        >>> [48.5, 38, 48.5, 38.5]
        
        Get pins for all the subswaths for two subswaths IW1 and IW2 processing:
        sbas.get_pins()
        >>> [47.5, 37.8, 47.5, 38.3, 48.5, 38, 48.5, 38.5]
        
        Get pins for all the subswaths for three subswaths IW1, IW2, IW3 processing:
        sbas.get_pins()
        >>> [47.5, 37.8, 47.5, 38.8, 48.5, 38, 48.5, 39, 49.5, 38.2, 49.5, 39.2]
```

```
def reframe_parallel(self, dates=None, n_jobs=-1)

		Reorder bursts from sequential scenes to cover the full orbit area between pins only.
		
    Args:
    		dates: optional list of dates. All the scenes processed when the argument is not defined.
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.

		Returns:
        None

    Examples:
    		sbas.reframe_parallel()
```

#### Build Interferogram (GMTSAR wrapper)

```
def stack_parallel(self, dates=None, n_jobs=-1):

		Stack and align scenes.
		
    Args:
    		dates: optional list of dates. All the scenes processed when the argument is not defined.
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.

		Returns:
        None

    Examples:
    		sbas.stack_parallel()
```

```
def intf_parallel(self, pairs, n_jobs=-1, wavelength=200, psize=32, func=None):

		Build interferograms for all the subswaths separately.
		
    Args:
    		pairs: list of dates pairs (baseline pairs).
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
    		wavelength: filtering wavelength, meters.
    		psize: patch size for modified Goldstein adaptive filter (power of two).
    		func: post-processing function usually used for decimation.

		Returns:
        None

    Examples:
    		For default 60m DEM resolution and other default parameters use command below:
    		pairs = [sbas.to_dataframe().index.unique()]
    		decimator = lambda dataarray: dataarray.coarsen({'y': 4, 'x': 4}, boundary='trim').mean()
    		sbas.intf_parallel(pairs, func=decimator)
```

#### Merging (GMTSAR wrapper)

```
def merge_parallel(self, pairs, grids = ['phasefilt', 'corr'], n_jobs=-1):

		Merge all the separate subswath interferograms into one large interferogram.
		
    Args:
    		pairs: list of dates pairs (baseline pairs).
    		grids: process one or the both interferogram related grids.
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
    		
		Returns:
        None

    Examples:
    		sbas.merge_parallel(pairs)
```

#### Manage landmask (PyGMTSAR original)

```
def set_landmask(self, landmask_filename):

    Define existing land mask file name for SBAS object.

    Args:
        landmask_filename: WGS84 NetCDF land mask file name.
        
    Returns:
        SBAS object.
        
    Examples:
    		sbas = sbas.set_landmask('data/landmask.nc')
    		Also, the same result is possible on SBAS initialisation:
    		sbas = SBAS(..., landmask_filename='data/landmask.nc')
```

```
def download_landmask(self, backend='GMT', debug=False):

    Download missed land mask file in geographic coordinates.

    Args:
        backend: argument to define processing backend. For now, only 'GMT' backend supported.
        debug: boolean flag to print debug information.
        
    Returns:
        None
        
    Examples:
    		sbas = sbas.download_landmask()
```

```
def get_landmask(self, geoloc=False, inverse_geocode=False, buffer_degrees=0.02):

		Return WGS84 land mask in geographic coordinates as Xarray Daraarray.
		
    Args:
    		geoloc: boolen flag to crop the land mask using the master scene approximate extent.
        inverse_: boolean flag to geocode land mask to radar coordinates. For this option the land mask
        cropped exactly to the interferogram grid.
        buffer_degrees: when geoloc=True add some buffer area around the master scene approximate extent.

    Returns:
        2D Xarray Daraarray object.

    Examples:
    		Get land mask in geographic coordinates:
        landmask_ll = sbas.get_landmask()
        
        Get land mask in radar coordinates:
        landmask_ra = sbas.get_landmask(inverse_geocode=True)
```

#### Unwrapping (SNAPHU wrapper)

```
def snaphu_config(self, defomax=0, **kwargs):

		Return SNAPHU text configuration to use it for unwrapping.
		
    Args:
    		defomax: define possible phase discontinuity [cycles].

    Returns:
        String.

    Examples:
    		Get default SHAPNU config:
        config = sbas.snaphu_config()
        
        Get custom SNAPHU config with added tiling options:
        config = sbas.snaphu_config(defomax=DEFOMAX, NTILEROW=1, NTILECOL=2, ROWOVRLP=200, COLOVRLP=200)
```

```
def unwrap_parallel(self, pairs, n_jobs=-1, threshold=None, conf=None, func=None, mask=None, conncomp=False):

		SNAPHU phase unwrapping.
		
    Args:
    		pairs: list of dates pairs (baseline pairs).
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
    		threshold: minimal valid correlation.
    		conf: SNAPHU text config.
    		func: post-processing function.
    		mask: Xarray Dataarray mask to exclude water sufraces (landmask) and other invalid areas.
    		conncomp: boolean flag to produce connected components grid.

		Returns:
        None

    Examples:
    		Simplest unwrapping:
    		sbas.unwrap_parallel(pairs)
    		
    		Filter low-coherence areas with common correlation threshold 0.75:
    		sbas.unwrap_parallel(pairs, threshold=0.075)
    		
    		Unwrap with coherence threshold 0.075 and fill NODATA gaps:
    		interpolator = lambda corr, unwrap: sbas.nearest_grid(unwrap).where(corr>=0)
    		sbas.unwrap_parallel(pairs, threshold=0.075, func=interpolator)
    		
    		Unwrap with coherence threshold 0.075 and apply land mask:
    		cleaner = lambda corr, unwrap: xr.where(corr>=0.075, unwrap, np.nan)
    		sbas.unwrap_parallel(pairs, threshold=0.075, mask=landmask_ra, func=cleaner)
    		
    		Unwrap with coherence threshold 0.075 and use SNAPHU tiling for faster processing and smaller RAM usage:
    		cleaner = lambda corr, unwrap: xr.where(corr>=0.075, unwrap, np.nan)
    		conf = sbas.PRM().snaphu_config(NTILEROW=1, NTILECOL=2, ROWOVRLP=200, COLOVRLP=200)
    		sbas.unwrap_parallel(pairs, n_jobs=1, threshold=0.075, func=cleaner, conf=conf)
```

#### Detrending (PyGMTSAR original)

```
def detrend(self, subswath, pair=None, wavelength=None, topo_ra=None, truncate=3.0, approximate=True, fit_intercept=True, fit_dem=True, fit_coords=True, debug=False):

		Detrend and get output for a single unwrapped interferogram combining topography and linear components removal plus Gaussian filtering.
		
    Args:
    		subswath: optional argument to define a subswath.
    		pairs: mandatory list of dates pairs (baseline pairs).
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
    		wavelength: optional cut-off wavelength for Gaussian filter in meters.
    		topo_ra: optional user-defined Xarray Dataarray topography grid instead of main one.
    		truncate: optional Gaussian filter window in sigmas.
    		approximate: optional boolean flag to use approximate Gaussian filter calculation for large gammas.
    		fit_intercept: optional boolean flag to drop mean value (plane).
    		fit_dem: optional boolean flag to detrend topography.
    		fit_coords: optional boolean flag to detrend linear coordinate components.
    		debug: boolean flag to print debug information.

		Returns:
        2D Xarray Dataarray

    Examples:
    		Simplest detrending:
    		unwrap_detrended = sbas.detrend(pair.values[0] if isinstance(pairs, pd.DataFrame) else pair[0])
    		Detrend ionospreric effects and solid Earth's tides on large area:
    		unwrap_detrended = sbas.detrend(pairs.values[0] wavelength=12000)
```

```
def detrend_parallel(self, pairs, n_jobs=-1, wavelength=None, topo_ra=None, truncate=3.0, approximate=True, fit_intercept=True, fit_dem=True, fit_coords=True):

		Detrend unwrapped interferograms combining topography and linear components removal plus Gaussian filtering. 
		
    Args:
    		pairs: list of dates pairs (baseline pairs).
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
    		wavelength: optional cut-off wavelength for Gaussian filter in meters.
    		topo_ra: optional user-defined Xarray Dataarray topography grid instead of main one.
    		truncate: optional Gaussian filter window in sigmas.
    		approximate: optional boolean flag to use approximate Gaussian filter calculation for large gammas.
    		fit_intercept: optional boolean flag to drop mean value (plane).
    		fit_dem: optional boolean flag to detrend topography.
    		fit_coords: optional boolean flag to detrend linear coordinate components.

		Returns:
        None

    Examples:
    		Simplest detrending:
    		sbas.detrend_parallel(pairs)
    		Detrend ionospreric effects and solid Earth's tides on large area:
    		sbas.detrend_parallel(pairs, wavelength=12000)
```

#### SBAS (GMTSAR wrapper)

```
def baseline_pairs(self, days=100, meters=150, invert=False):

		Generate SBAS baseline pairs. 
		
    Args:
    		days: maximum baseline in days.
    		meters: maximum perpendicular baseline in meters.
    		invert: boolean flag to invert all the pairs.

		Returns:
        Pandas Dataframe.

    Examples:
    		sbas.baseline_pairs()
```

```
def sbas(self, pairs, smooth=0, atm=0, debug=False):

		SBAS unwrapped interferograms timeseries analysis, see GMTSAR documentation for the details.
		
    Args:
    		pairs: list of dates pairs (baseline pairs).
    		smooth: optional smoothing factor.
    		atm: optional atmospheric correction factor.
    		debug: boolean flag to print debug information.

		Returns:
        None

    Examples:
    		sbas.sbas()
```



#### SBAS (PyGMTSAR original)

```
def sbas_parallel(self, pairs, mask=None, detrended=True, n_jobs=-1):

		SBAS unwrapped and detrended interferograms timeseries analysis using pixel-wise correlation-weighted least squares approach. 
		
    Args:
    		pairs: list of dates pairs (baseline pairs).
    		mask: Xarray Dataarray mask for valid pixels.
    		detrended: optional boolean flag to use detrended unwrapped phase or unwrapped phase itself.
    		n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.

		Returns:
        None

    Examples:
    		sbas.sbas_parallel()
```

#### Incidence angle (PyGMTSAR original)

```
def incidence_angle(self, subswath=None):

		Compute incidence angle grid in geographic coordinates. 
		
    Args:
    		subswath: optional subswath number.

		Returns:
        2D Xarray Dataaarray.

    Examples:
    		da_ll = sbas.incidence_angle()
```

#### Displacements (PyGMTSAR original)

```
def los_displacement_mm(self, unwraps):

		Compute LOS displacement in millimeters. 
		
    Args:
    		unwraps: unwrapped phase 2D grid or 3D grids stack.

		Returns:
        2D or 3D Xarray Dataaarray.

    Examples:
    		Calculate LOS displacement for unwrapped phase grids in radar coordinates:
    		unwraps_ra = sbas.open_grids(pairs, 'unwrap')
    		los_disp_ra = sbas.los_displacement_mm(unwraps_ra)
    		or the same code in one line
    		los_disp_ra = sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm)
    		Note: here "func" argument for open_grids() function reduces the code to a single command.
    		
    		Calculate LOS displacement for detrended unwrapped phase grids in geographic coordinates:
    		detrend_ll = sbas.open_grids(pairs, 'detrend', geocode=True)
    		los_disp_ll = sbas.los_displacement_mm(detrend_ll)
    		or the same code in one line
    		los_disp_ll = sbas.open_grids(pairs, 'detrend', geocode=True, func=sbas.los_displacement_mm)
    		Note: here "func" argument for open_grids() function reduces the code to a single command.
```

```
def vertical_displacement_mm(self, unwraps):

		Compute vertical displacement in millimeters in geographic coordinates. 
		
    Args:
    		unwraps: unwrapped phase 2D grid or 3D grids stack.

		Returns:
        2D or 3D Xarray Dataaarray.

    Examples:
    		Calculate vertical displacement for unwrapped phase grids in geographic coordinates:
    		unwraps_ll = sbas.open_grids(pairs, 'unwrap', geocode=True)
        vert_disp_mm = sbas.los_displacement_mm(unwraps_ll)
        
        Calculate vertical displacement for detrended unwrapped phase grids in geographic coordinates:
    		vert_disp_mm = sbas.open_grids(pairs, 'detrend', geocode=True, func=sbas.vertical_displacement_mm)
    		Note: here "func" argument for open_grids() function reduces the code to a single command.
```

```
def eastwest_displacement_mm(self, unwraps):

		Compute East-West displacement in millimeters. 
		
    Args:
    		unwraps: unwrapped phase 2D grid or 3D grids stack.

		Returns:
        2D or 3D Xarray Dataaarray.

    Examples:
    		Calculate East-West displacement for unwrapped phase grids in geographic coordinates:
    		unwraps_ll = sbas.open_grids(pairs, 'unwrap', geocode=True)
        ew_disp_mm = sbas.eastwest_displacement_mm(unwraps_ll)
        
    		Calculate East-West displacement for detrended unwrapped phase grids in geographic coordinates:
    		ew_disp_mm = sbas.open_grids(pairs, 'detrend', geocode=True, func=sbas.eastwest_displacement_mm)
    		Note: here "func" argument for open_grids() function reduces the code to a single command.
```

#### Data output (PyGMTSAR original)

```
def open_grids(self, pairs, name, geocode=False, mask=None, func=None,
               crop_valid=False, add_subswath=True,  n_jobs=-1):

    Lazy open PyGMTSAR produced NetCDF grids as 3D data cube and apply a set of functions on it.

    Args:
        pairs: list of dates pairs (baseline pairs) or dates.
        name: grid name one of the list: 'phasefilt', 'corr', 'unwrap', 'detrend', 'disp'.
        geocode: optional boolean flag to geocode the grid to geographic coordinates.
        mask: optionalXarray Dataarray mask to exclude invalid areas.
        func: optional function name or a list of function names to apply too the each one 2D input grid,
        crop_valid: optional boolean flag to crop valid grid extent only.
        add_subswath: optional boolean flag to add subswath name as 'Fn' to grid names.        
        n_jobs: number of parallel processing jobs. n_jobs=-1 means all the processor cores used.

    Returns:
        3D Xarray Dataarray Dask Dataaarray.

    Examples:
        Open unwrapped phase grids in radar coordinates:
        unwraps_ra = sbas.open_grids(pairs, 'unwrap')

        Open PyGMTSAR SBAS grids in radar coordinates and calculate LOS displacement and fill NODATA areas:
        dates = np.unique(pairs.values.flatten() if isinstance(pairs, pd.DataFrame) else pairs.flatten())
        disps = sbas.open_grids(dates, 'disp', func=[sbas.los_displacement_mm, sbas.nearest_grid])

        Open GMTSAR SBAS grids in radar coordinates using compatibility option add_subswath=False:
        sbas.open_grids(sbas.df.index, 'disp', func=sbas.nearest_grid, add_subswath=False)

        Calculate LOS displacement for unwrapped phase grids in radar coordinates:
        unwraps_ra = sbas.open_grids(pairs, 'unwrap')
        los_disp_ra = sbas.los_displacement_mm(unwraps_ra)
        or the same code in one line
        los_disp_ra = sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm)
        Note: here "func" argument for open_grids() function reduces the code to a single command.

        Calculate LOS displacement for detrended unwrapped phase grids in geographic coordinates:
        detrend_ll = sbas.open_grids(pairs, 'detrend', geocode=True)
        los_disp_ll = sbas.los_displacement_mm(detrend_ll)
        or the same code in one line
        los_disp_ll = sbas.open_grids(pairs, 'detrend', geocode=True, func=sbas.los_displacement_mm)
        Note: here "func" argument for open_grids() function reduces the code to a single command.

        Open GMTSAR SBAS velocity grid in radar coordinates:
        vel = sbas.open_grids(None, 'vel', add_subswath=False)
        Open GMTSAR SBAS velocity grid in geographic coordinates:
        vel = sbas.open_grids(None, 'vel', geocode=True, add_subswath=False)
```

#### Geocoding (PyGMTSAR original)

```
def intf_ra2ll(self, subswath=None, grids=None):
    
		Geocoding function based on interferogram geocode matrix.
		
	Args:
			subswath: optional argument to define a subswath.
			grids: mandatory 2D or 3D Xarray Dataaarray in radar coordinates.      
	
		Returns:
	    2D or 3D Xarray Dataaarray.
	
	Examples:
			Geocode 3D unwrapped phase grids stack:
	    unwraps_ll = sbas.intf_ra2ll(sbas.open_grids(pairs, 'unwrap'))
	    or use "geocode" option instead:
	    unwraps_ll = sbas.open_grids(pairs, 'unwrap', geocode=True)
```

```
def intf_ll2ra(self, subswath=None, grids=None):
    
		Inverse geocoding function based on interferogram inverse geocode matrix.
		
	Args:
			subswath: optional argument to define a subswath.
			grids: mandatory 2D or 3D Xarray Dataaarray in geographic coordinates.      
	
		Returns:
	    2D or 3D Xarray Dataaarray.
	
	Examples:
			Inverse geocode 3D unwrapped phase grids stack:
			unwraps_ll = sbas.open_grids(pairs, 'unwrap', geocode=True)
	    unwraps = sbas.intf_ll2ra(unwraps_ll)
```

#### Backup and restore (PyGMTSAR original)
```
def dump(self, to_path=None):
    
		Dump SBAS object state to pickle file (SBAS.pickle in the processing directory by default).
		
	Args:
			to_path: optional path to output dump file.      
	
		Returns:
	    None
	
	Examples:
			Dump current state to default dump file in the processing directory:
			sbas.dump()
```

```
def restore(from_path=None):
    
		Restore SBAS object state from pickle file (SBAS.pickle in the processing directory by default).
		
	Args:
			from_path: optional path to input dump file.      
	
		Returns:
	    None
	
	Examples:
			Restore current state from default dump file in the processing directory:
			sbas.restore()
```

```
def backup(self, backup_dir):
    
		Backup framed SBAS scenes, orbits, DEM and landmask files to build a minimal reproducible dataset.
		
	Args:
			backup_dir: backup directory.      
	
		Returns:
	    None
	
	Examples:
			Backup to the specified directory:
			sbas.backup('backup')
			Open the backup for the reproducible run defining it as new data directory:
			sbas = SBAS('backup', 'backup/DEM_WGS84.nc', 'raw')
```

#### Helpers (PyGMTSAR original)

```
def pixel_spacing(self, grid=(1, 4)):
    
		Compute pixel size in meters for the default processing grid or a current one.
		
	Args:
			grid: a pair of x,y grid decimation coefficients or 2D or 3D Xarray Dataarray.      
	
		Returns:
	    Pair of float numbers.
	
	Examples:
			Get default pixel spacing:
			sbas.pixel_spacing()
			>>> (14.0, 14.8)
			
			Get unwrapped phase grid actual pixel spacing for interferogram decimation {'y': 2, 'x': 2}:
			sbas.pixel_spacing(unwraps)
			>>> (27.9, 29.5)
```

```
def nearest_grid(in_grid, search_radius_pixels=300):
    
		Nearest neighbor interpolation.
		
	Args:
			in_grid: 2D Xarray Dataarray. 
      search_radius_pixels: optional interpolation distance in pixels.
	
		Returns:
	    2D Xarray Dataarray.
	
	Examples:
			Fill NODATA gaps in the specified grid:
			sbas.search_radius_pixels(da)
```

```
def cropna(self, das):
    
		Crop raster valid extent removing all rows and columns containing NODATA values only.
		
	Args:
			das: 2D or 3D Xarray Dataarray. 
	
		Returns:
	    2D or 3D Xarray Dataarray.
	
	Examples:
			Crop valid raster extent only:
			sbas.cropna(unwraps)
```

```
def as_geo(self, da):
    
		Add geospatial attributes (CRS and spatial dimentions) to allow RioXarray raster operations.
		
	Args:
			da: 2D or 3D Xarray Dataarray. 
	
		Returns:
	    2D or 3D Xarray Dataarray.
	
	Examples:
			Convert raster to geospatial and mask it by a Shapely vector geometry:
			sbas.as_geo(unwraps).rio.clip([geometry])
```


## Installation

PyGMTSAR project includes GMTSAR binary tools plus Python library and the installation requires as the binaries as the library. Miss the binary installation in case when you already have the recent GMTSAR installed.

### Install PyGMTSAR Python library on MacOS (Apple Silicon and Intel)

[Homebrew](https://brew.sh) package manager is required and it should be installed first, follow the link to do it. As the installation step Apple command line developer tools needt to be installed:

```
xcode-select --install
```

The commands below tested on Intel BigSur and Apple Silicon Monterey systems and perhaps work on older MacOS versions too. Python 3.10 is strictly required for parallel computations on Apple Silicon and it's recommended on Intel Macs too while thta's possible to use Python 3.7+ as on Ubuntu 18.04 (bionic). 

```
# Python 3.10 itself
brew install python@3.10
# NetCDF supporting packages
brew install hdf5 netcdf
# install system GIS packages
brew install proj gdal
```

NumPy library can be significantly accelerated (~4x) using [Apple Accelerate Framework](https://developer.apple.com/documentation/accelerate) for using custom compilation as

```
python3 -m pip uninstall -y numpy
python3 -m pip install cython pybind11
# numba requires recent numpy 1.22.4
python3 -m pip install --no-binary :all: --no-use-pep517 numpy==1.22.4
python3 -m pip install numba
# check the output for string '-Wl,Accelerate' to be sure the acceleration works well
python3 -c "import numpy; numpy.show_config()"
```

And install some exact library versions to prevent dependencies conflicts:

```
# this version requires proj (8+) including proj.h
python3 -m pip install cartopy==0.20.2
python3 -m pip install xarray==0.19.0
python3 -m pip install scipy==1.8.1
```

All other packages installation is straightforward:

```
python3 -m pip install \
    h5py netcdf4 h5netcdf \
    rasterio scikit-image \
    distributed zarr nc-time-axis \
    matplotlib seaborn geoviews hvplot datashader bokeh
```

Note: SciPy does not support [Apple Accelerate Framework](https://developer.apple.com/documentation/accelerate) on MacOS and SciPy custom compilation is useless.

Install the library using PIP packages manager as

```
python3 -m pip install pygmtsar
```

Optionally check the installation from you Python shell or Jupyter notebook:

```
from pygmtsar import SBAS
```

### Install PyGMTSAR Python library on Debian 10

Install PyGMTSAR dependencies including system packages and Python libraries:

```
apt install -y python3-pip
python3 -m pip install sentineleof elevation
eio selfcheck

# for python tools and sbas command line arguments calculation
apt install -y gdal-bin python-gdal python3-netcdf4 python3-scipy bc
# for additional Python-coded utilities
python3 -m pip install xarray numpy scipy pytest --upgrade
# for notebooks
apt -y install libgmt5 netcdf-bin python3-pip python3-netcdf4 python3-scipy python3-dev
apt -y install gdal-bin python3-gdal libcharls2 libgdal-dev libproj-dev proj-data proj-bin libgeos-dev
```

Install the library using PIP packages manager as

```
python3 -m pip install pygmtsar
```

Optionally check the installation from you Python shell or Jupyter notebook:

```
from pygmtsar import SBAS
```

### Install PyGMTSAR Python library on Ubuntu 18.04 (bionic)

Install PyGMTSAR dependencies including system packages and Python libraries:

```
# use the old versions vbecause the modern ones does not work properly on the outdated Ubuntu
python3 -m pip install cartopy==0.19.0.post1  1>/dev/null 2>/dev/null
python3 -m pip install install xarray==0.19.0 1>/dev/null 2>/dev/null
python3 -m pip install install scipy==1.7.1   1>/dev/null 2>/dev/null

python3 -m pip install \
    h5py netcdf4 h5netcdf \
    rasterio scikit-image \
    distributed zarr nc-time-axis \
    matplotlib seaborn geoviews hvplot datashader bokeh
```

Install the library using PIP packages manager as

```
pip3 install pygmtsar
```

Optionally check the installation from you Python shell or Jupyter notebook:

```
from pygmtsar import SBAS
```

### Install PyGMTSAR Python library on some other operation system

Install the library using PIP packages manager as

```
python3 -m pip install pygmtsar
```

This command will try to install all the required Python libraries but some system packages might be required to install manually. Probably you know how to do it when you use a rare operation system.

Optionally check the installation from you Python shell or Jupyter notebook:

```
from pygmtsar import SBAS
```

### Install GMTSAR binaries on MacOS (Apple Silicon and Intel)

[Homebrew](https://brew.sh) package manager is required and it should be installed first, follow the link to do it. As the installation step Apple command line developer tools needt to be installed:

```
xcode-select --install
```

And some system dependencies required:

```
brew install wget libtiff hdf5 gmt ghostscript autoconf
```

After that the installation is straightforward, run the same commands for BigSur on Intel chips and Monterey on Apple Silicon chips:

```
# create installation directory
sudo mkdir /usr/local/GMTSAR
sudo chown $(whoami) /usr/local/GMTSAR
# install recent GMTSAR
cd /usr/local
git clone --branch master https://github.com/gmtsar/gmtsar GMTSAR
cd /usr/local/GMTSAR
autoconf
./configure --with-orbits-dir=/tmp
make
make install
```

Note: that's possible to install GMTSAR to /opt directory instead on Apple Silicon Macs.

### Install GMTSAR binaries on Debian 10

There is the cloud initialisation scripit [GMTSAR.install.debian10.sh](https://github.com/mobigroup/gmtsar/blob/master/gmtsar/sh/GMTSAR.install.debian10.sh) to install and configure all the dependencies and GMTSAR on cloud Debian 10 hosts.

```
GMTSAR=/usr/local/GMTSAR
GIT=https://github.com/gmtsar/gmtsar
BRANCH=master

# prepare system
apt update
apt -y upgrade
apt install -y locales
# Uncomment en_US.UTF-8 for inclusion in generation
sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen
# Generate locale
locale-gen en_US.UTF-8
# Export env vars
echo "export LC_ALL=en_US.UTF-8" >> /root/.bashrc
echo "export LANG=en_US.UTF-8" >> /root/.bashrc
echo "export LANGUAGE=en_US.UTF-8" >> /root/.bashrc

# https://github.com/gmtsar/gmtsar/wiki/GMTSAR-Wiki-Page
apt install -y csh subversion autoconf libtiff5-dev libhdf5-dev wget
apt install -y liblapack-dev
apt install -y gfortran
apt install -y g++
apt install -y libgmt-dev
apt install -y gmt-dcw gmt-gshhg
# gmt-gshhg-full should be installed automatically (it is required to use GMTSAR landmask)
apt install -y gmt
# fix for missed "gs" utility and git and make tools
apt install -y git make

cd $(dirname "$GMTSAR")
git clone --branch "$BRANCH" "$GIT" GMTSAR
cd GMTSAR
autoconf
./configure --with-orbits-dir=/tmp
make
make install
```

### Install GMTSAR binaries on Ubuntu 18.04 (bionic)

This operation system is installed on Google Colab platform and the Live Google Colab notebooks include the installation commands.

```
apt install -y csh autoconf gfortran \
    libtiff5-dev libhdf5-dev liblapack-dev libgmt-dev gmt-dcw gmt-gshhg gmt > /dev/null
cd /usr/local && git clone --branch master https://github.com/gmtsar/gmtsar GMTSAR > /dev/null
cd /usr/local/GMTSAR && autoconf > /dev/null
cd /usr/local/GMTSAR && ./configure --with-orbits-dir=/tmp > /dev/null
cd /usr/local/GMTSAR && make 1>/dev/null 2>/dev/null
cd /usr/local/GMTSAR && make install >/dev/null
```

### Install GMTSAR binaries on some other operation system

GMTSAR compilation is possible on any POSIX-compatible Unix/Linux/BSD system while some system packages should be installed. See for more recipes [GMTSAR Wiki Page](https://github.com/gmtsar/gmtsar/wiki/GMTSAR-Wiki-Page)



### @ Alexey Pechnikov, September, 2022

<!-- Setup title and logo -->
<script>
document.querySelector("header h1").textContent = 'PyGMTSAR'
document.querySelector("header p").textContent = 'Easy and Fast Satellite Interferometry For Everyone'
//document.querySelector("header p a").textContent = 'View on GitHub mobigroup/gmtsar'
this.img = document.createElement("img");
this.img.style.cssText = 'width: 200px';
this.img.src = "https://user-images.githubusercontent.com/7342379/190416030-5388bf04-e322-4616-9d33-adf84247d976.png";
src = document.querySelector("p.view");
src.appendChild(this.img);
</script>

<style>
div.wrapper {width: 960px;}
section {width: 720px;}
header {width: 230px;}
footer {width: 230px;}
</style>
