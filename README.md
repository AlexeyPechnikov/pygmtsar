__INSTRUCTIONS FOR INSTALLING AND RUNNING GMTSAR__
----------------------------------------------

 David T. Sandwell  -  dsandwell@ucsd.edu
 
 Eric Xu            -  xix016@ucsd.edu
 
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
February 13, 2010 Code available as simple tar file;
April     4, 2013 Converted to use GMT5 API, available with SVN;
July     28, 2019 Available with github

__INSTALL__

1) Go to the GMTSAR WIKI and follow the instructions.
       http://gmt.soest.hawaii.edu/projects/gmt5sar/wiki

2) Download orbit files for just ERS and ENVISAT. The other satellites have orbits distributed with the data.
   Put it anywhere you like; let's refer to this dir as <orbitsdir>
       https://topex.ucsd.edu/gmtsar/downloads/

3) Go to the gmtsar directory and enter:

       autoconf

4) Configure GMTSAR:	
   Run the configure script.  To see all options, run
   
       ./configure --help
	
   Most users will simply run
   
       ./configure --with-orbits-dir=<orbitsdir> --prefix=<installdir>
	
   where <orbitsdir> was defined in step (3) and <installdir> is where you want to place all the executables.
   For example, you might run
	
       ./configure --with-orbits-dir=/usr/local/orbits

5) To build all executables, type
       make
   Assuming that went well, you can install executables, scripts, and shared data by running
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
