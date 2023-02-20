[![GMTSAR tests](https://github.com/gmtsar/gmtsar/actions/workflows/gmtsar.yml/badge.svg)](https://github.com/gmtsar/gmtsar/actions/workflows/gmtsar.yml)

__INSTRUCTIONS FOR INSTALLING AND RUNNING GMTSAR__
----------------------------------------------

__INSTALL__

1) Go to the GMTSAR WIKI and follow the instructions to install dependencies and download the GMTSAR package.
       https://github.com/gmtsar/gmtsar/wiki/GMTSAR-Wiki-Page

2) Download orbit files for just ERS and ENVISAT. The other satellites have orbits distributed with the data.
       https://topex.ucsd.edu/gmtsar/downloads/

3) Go to the gmtsar directory and enter:

       autoconf

4) Configure GMTSAR:	
   Run the configure script.  To see all options, run
   
       ./configure --help
	
   Most users will simply run
   
       ./configure --with-orbits-dir=<orbitsdir> --prefix=<installdir>
	
   For example, you might run
	
       ./configure --with-orbits-dir=/usr/local/orbits

5) To build all executables, type

       make
       make install

6) test the commands: gmt, esarp, xcorr, conv, gmtsar.csh, etc.
   If using C-shell you may have to type rehash first. 
   If this does not work then make sure the <installdir> is in your system $PATH or $path.

__RUN__

1) GET DATA. There is an example data set at our website:
       https://topex.ucsd.edu/gmtsar/downloads/

Uncompress the file and then unpack with tar.

2) ORGANIZE the DISK. The standard GMTSAR run has the following directories. 

  raw    SLC    topo   intf

raw - contains the original data.

SLC - contains the single look complex images derived from the raw data.

topo - contains a digital elevation model for the area in geographic and ultimately radar coordinates.

intf - contains subdirectories with the possible interferograms (only 1 for this example).

3) MAKE DEM 

Go to the following web site and construct a file dem.grd that encloses the SAR frame.
      https://topex.ucsd.edu/gmtsar/demgen

Place the file dem.grd in the /topo directory. 

4) EXAMPLE RUNS

There are complete examples for each data type at the following web site:
      https://topex.ucsd.edu/gmtsar/downloads/

To process a pair interferogram, try `p2p_processing.csh`, or `p2p_S1_TOPS_Frame.csh` for Sentinel-1 TOPS data and `p2p_ALOS2_SCAN_Frame.csh` for ALOS-2 ScanSAR data. For batch processing, see instructions in the GMTSAR documentation (http://topex.ucsd.edu/gmtsar/tar/GMTSAR_2ND_TEX.pdf) or the Sentinel-1 time-series recipe (http://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_5.pdf).
