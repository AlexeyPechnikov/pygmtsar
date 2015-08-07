#!/bin/csh -f
#       $Id$
# Xiaopeng Tong Nov 23 2009
#
unset noclobber
# 
# Converts a complex SLC file to a real amplitude grd
# file using optional filter and a PRM file
#
# define the filters
#
  set sharedir = `gmtsar_sharedir.csh`
  set fil1 = $sharedir/filters/gauss5x3
  set fil2 = $sharedir/filters/gauss9x5 
#
# check for number of arguments
#
  if ($#argv != 3 ) then
    echo ""
    echo "Usage: slc2amp.csh filein.PRM rng_dec fileout.grd " 
    echo ""
    echo "       rng_dec is range decimation"
    echo "       e.g. 1 for ERS ENVISAT ALOS FBD" 
    echo "            2 for ALOS FBS " 
    echo "            4 for TSX"
    echo ""
    echo "Example: slc2amp.csh IMG-HH-ALPSRP055750660-H1.0__A.PRM 2 amp-ALPSRP055750660-H1.0__A.grd"
    echo ""
    exit 1
  endif 
#
# filter the amplitudes done in conv
# check the input and output filename before 
#
  if ((($1 =~ *PRM*) || ($1 =~ *prm*)) && ($3 =~ *grd*)) then
    echo " range decimation is:" $2
    conv 4 $2 $fil1 $1 $3
  else 
    echo "slc2amp.csh"
    echo "wrong filename" 
    exit 1
  endif
# 
# get the zmin and zmax value
#
  gmt grdmath $3 1 MUL = $3
