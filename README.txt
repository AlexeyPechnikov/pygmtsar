*****************************************************************************************************
#	$Id: README.TXT 211 2015-08-07 15:49:42Z pwessel $
#
INSTRUCTIONS FOR INSTALLING AND RUNNING GMT5SAR 

Copyright (c) 2009-2015
David T. Sandwell  -  dsandwell@ucsd.edu
Rob Mellors        -  rmellors@geology.sdsu.edu
Xiaopeng Tong      -  xitong@ucsd.edu
Meng Wei           -  mwei@ucsd.edu
Paul Wessel        -  pwessel@hawaii.edu

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 3 of the LICENSE.TXT

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

Significant Modifications:
February 13, 2010
March    11, 2010
March    25, 2010
June      8, 2010
September 27, 2010
April 4, 2013	Converted to use GMT5 and redone the entire install procedure.  Linking with GMT5
	means GMT5SAR will use the FFTs recognized by GMT5: Accelerate Framework on OS X, FFTW,
	KISS FFT, or Brenner Fortran-translated FFT. GMT5 selects the fastest FFT given the
	dimensions and the availability of libraries.

*****************************************************************************************************

1) Obtain GMT5 from SOEST and install [see "Obtaining GMT" on gmt.soest.hawaii.edu].
   Make sure all the required libraries (e.g., netcdf, GDAL, PCRE) have been installed
   before you start the GMT build/install procedure.

2) Go to the GMT5SAR WIKI and follow the instructions.
   http://gmt.soest.hawaii.edu/projects/gmt5sar/wiki

3) Download orbit files for ERS and ENVISAT. No need for ALOS orbit data because it is built in raw file.
   The file is > 1G and contains all the available orbit for ERS and ENVISAT. It takes a while to download
   but you only need to download it once. Put it anywhere you like; let's refer to this dir as <orbitsdir>

   http://topex.ucsd.edu/gmtsar/downloads

4) Obtain GMT5SAR subversion. Go to the GMT5SAR WIKI and follow the instructions.
   http://gmt.soest.hawaii.edu/projects/gmt5sar/wiki
     
	autoconf

   (When we install you can specify a system directory such as /usr/local, assuming you have permission).

4) Configure GMT5SAR:	
   Run the configure script.  To see all options, run

	./configure --help
	
   Most users will simply run

	./configure --with-orbits-dir=<orbitsdir> --prefix=<installdir>
	
   where <orbitsdir> was defined in step (3) and <installdir> is where you want to place all the executables.
   For example, you might run

	./configure --with-orbits-dir=/usr/local/orbits

   If you are a developer you may also want to add --enable-debug so you can
   run the programs in a visual debugger.
   Note: If you ever decide to move the orbits directory then you must also reconfigure GMT5SAR.
   
5) To build all executables, type
	make
   Assuming that went well, you can install executables, scripts, and shared data by running
	make install
   or
	sudo make install
   if you are writing to a system directory that requires admin privileges.
   After this step you can remove your staging directory/tar-files if you like.

6) test the commands: gmt, esarp, xcorr, conv, gmtsar.csh, etc.
   If using C-shell you may have to type rehash first. 
   If this does not work then make sure the <installdir> is in your system $PATH or $path.

*****************************************************************************************************
RUN

1) GET DATA. There is an example data set at our website:

http://topex.ucsd.edu/gmtsar/downloads/

Uncompress the file and then unpack with tar.

2) ORGANIZE the DISK. The standard GMTSAR run has the following directories. 

  raw 	SLC 	topo   intf

raw - contains the original data.  

SLC - contains the single look complex images derived from the raw data.

topo - contains a digital elevation model for the area in geographic and ultimately radar coordinates.

intf - contains subdirectories with the possible interferograms (only 1 for this example).

In raw, original data have special name format. If you got data that have different name system, you need to change the name according to the following:
 
For ALOS we just need the IMG- and LED-files.  Here is an example.
-rw-r--r--  1 sandwell  14  382547520 Aug 28 15:29 IMG-HH-ALPSRP227730640-H1.0__A
-rw-r--r--  1 sandwell  14   12506972 Aug 28 15:28 LED-ALPSRP227730640-H1.0__A
-rw-r--r--  1 sandwell  14   12506972 Aug 28 15:28 LED-ALPSRP207600640-H1.0__A
-rw-r--r--  1 sandwell  14  747383820 Aug 28 15:28 IMG-HH-ALPSRP207600640-H1.0__A

Make sure the images are all from the same orbit direction, same track, and same frame.  The filename 
ALPSRP207600640-H1.0__A  can be decomposed as:

__A - ascending
0640 - frame number
20760 - orbit number   To check that the data are from the same track, the difference between orbit numbers should be divisible by 671.   In this case there are 7,  46-day cycles between the acquisitions.

For ENVISAT we just need .baq file. Here is an example.
ENV1_2_084_2943_42222.baq
ENV1_2_084_2943_42723.baq

The filename means:
ENVISAT _ 2 _ track number _ frame number _ orbit number

For ERS we just need .dat and .ldr file. Here is an example.
e2_127_2907_23027.ldr
e2_127_2907_23027.dat
e2_127_2907_23528.ldr
e2_127_2907_23528.dat

The filename means:
ERS2_track number_frame number_orbit number.
.dat - raw data file
.ldr - leader file

3) MAKE DEM 

Go to the following web site and construct a file dem.grd that encloses the SAR frame.
http://topex.ucsd.edu/gmtsar/demgen

Place the file dem.grd in the /topo directory. 

4) DO EVERYTHING

For ALOS:
p2p_ALOS.csh IMG-HH-ALPSRP207600640-H1.0__A IMG-HH-ALPSRP227730640-H1.0__A configure.txt

The file configure.txt should be edited to set a number of parameters including the starting point for t
he InSAR processing.

For ENVISAT:
p2p_ENVI.csh ENV1_2_084_2943_42222 ENV1_2_084_2943_42723 dem.grd

For ERS:
p2p_ERS.csh e2_127_2907_23027 e2_127_2907_23528 dem.grd

*********************************************************************************************************

RUN - STEP-BY-STEP INSTRUCTIONS

Taking ALOS for example but applicable to ERS and ENVISAT.

1) PREPROCESS the raw data

    cd raw
    ls IMG* >> data.in

Edit the data.in file and place the master in the first line. 
Also one can set the -radius and -near_range on the first line to have this frame match other frames along the same track.
When this command is done you will have PRM and raw files for every scene in the list data.in, all in the same format and
geometry. Execute the command.

   pre_proc_batch.csh ALOS data.in

2) ALIGN the SLC images.  Link the raw data into the SLC directory.
    cd SLC
    cp ../raw/*.PRM .
    ln -s ../raw/*.raw .
    ln -s ../raw/LED* .

Align the images
   
    align.csh ALOS IMG-HH-ALPSRP207600640-H1.0__A IMG-HH-ALPSRP227730640-H1.0__A 

This will keep your computer busy for a while.  The first image is called the master and the second is a slave.  One can
align a slave to the master and then use that slave as a surrogate master.  This is useful for alignment of a large
stack of data having a large spread in the baseline versus time plot.  When this is done you will have 2 output files
for each scene, a PRM-file and a Single Look Complex (SLC-file).   The PRM-files are ascii text files containing the
parameters needed for InSAR processing.

master
IMG-HH-ALPSRP207600640-H1.0__A.PRM
IMG-HH-ALPSRP207600640-H1.0__A.SLC

slave
IMG-HH-ALPSRP227730640-H1.0__A.PRM
IMG-HH-ALPSRP227730640-H1.0__A.SLC

(Look inside the script align.csh to see what it does.)

3) MAKE the topo_ra.grd

   cd topo
   cp ../SLC/IMG-HH-ALPSRP207600640-H1.0__A.PRM master.PRM
   ln -s ../raw/LED-ALPSRP207600640-H1.0__A .

Next you will need  a file dem.grd.  Thus can be created at:
http://topex.ucsd.edu/gmtsar/demgen

Once you have the dem.grd in lon/lat coordinates you can convert it to radar coordinates using the following command. 

  dem2topo_ra.csh master.PRM dem.grd

This creates a file called topo_ra.grd.  The script also creates postscript images of the dem.ps and topo_ra.ps that can be viewed.
The script also creates a file trans.dat that has the complete mapping from lon,lat,topo to range,azimuth.  This same file will
be used later for geocoding, (i.e., converting the range, azimuth grids back into lon, lat space).

4) INTERFEROGRAM

  cd intf
  mkdir 20760_22773  

Note these are the orbit numbers of the reference and repeat images.  One could also use dates for the directory name. 

   cd 20760_22773
   ln -s ../../raw/LED-ALPSRP207600640-H1.0__A .
   ln -s ../../raw/LED-ALPSRP227730640-H1.0__A .
   cp ../../SLC/IMG-HH-ALPSRP207600640-H1.0__A.PRM .
   cp ../../SLC/IMG-HH-ALPSRP227730640-H1.0__A.PRM .
   ln -s ../../SLC/IMG-HH-ALPSRP207600640-H1.0__A.SLC .
   ln -s ../../SLC/IMG-HH-ALPSRP227730640-H1.0__A.SLC .
   ln -s ../../topo/topo_ra.grd .

All these files are needed to be linked or copied to make the interferogram.  Of course this could all be done with a script.
Now make the interferogram.  While it is running look inside the script intf.csh to see the grdmath.  Also, if you have any
kind of error, delete everything and start over.

  intf.csh IMG-HH-ALPSRP207600640-H1.0__A.PRM IMG-HH-ALPSRP227730640-H1.0__A.PRM  -topo topo_ra.grd

The output will be three grd  and postscript files.

display_amp.grd, display_amp.ps  - amplitude of interferogram
phase.grd, phase.ps              - phase of interferogram
corr.grd, corr.ps                - correlation of interferogram

The phase measures the displacement of the repeat images relative to the reference images.
The "reference and repeat" is a separate definition from the "master and slave" (refer to the document).
Note that the phase is relative measurement so it's not important whether the pixel values are negative or positive.
Phase increase means that the ground is moving away from the radar (range increasing);
phase decrease means that the ground is moving toward the radar (range decreasing).

5) FILTER INTERFEROGRAM

   filter.csh IMG-HH-ALPSRP207600640-H1.0__A.PRM IMG-HH-ALPSRP227730640-H1.0__A.PRM gauss_alos_200m 2

 The gaussian filter and decimation for the amplitude and phase images are the same.
 filter is the  name of the filter.
 decimation control the size of the amplitude and phase images. It is either 1 or 2.
 Set the decimation to be 1 if you want higher resolution images.
 Set the decimation to be 2 if you want images with smaller file size.

The output will be masked and filtered phase as well as phase gradients (xphase.grd and yphase.grd).

6) UNWRAP PHASE

   cd intf/20760_22773
   shaphu.csh .10  

7) GEOCODE

   cd intf/20760_22773
   ln -s ../../topo/trans.dat .
   geocode.csh .10

This command will make postscript and KML images of the phase, correlation, and display amplitude.
The argument 0.15 is used to mask the phase when the coherence is less than 0.15.  Use your GMT skills to
improve on these maps.  One can combine the phase (color) and amplitude (shading) in gmt grdimage or combine
the phase (color) and gmt grdgradient dem.grd (shading) to make interesting plots.

8) SHIFT THE TOPOPHASE [OPTIONAL]

The topo_ra.grd may not be perfectly aligned with the master SLC. This can be corrected by shifting
the topo_ra by 1 or 2 pixels, usually in the range direction.  

   cd SLC
   slc2amp.csh IMG-HH-ALPSRP207600640-H1.0__A.PRM amp_master.grd
   cd ../topo
   ln -s ../SLC/amp_master.grd .
   offset_topo amp_master.grd topo_ra.grd 0 0 5 topo_shift.grd

The last line of the output from this program shows the shift needed to maximize the cross correlation
between the amplitude of the master image and the range gradient of the topo_ra.grd.  The output file
is the topography shifted by an integer number of pixels.  Now go back to step 4) and use the topo_ra_shift.grd
instead of the topo_ra.grd to remake the interferogram.
