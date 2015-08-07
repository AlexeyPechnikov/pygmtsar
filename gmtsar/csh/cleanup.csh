#!/bin/csh -f
#       $Id$
#
#  David T. Sandwell, MAR 11, 2010
#
# Clean the disk area in preparation for process2pass.csh
# This should be run in the top directory.  An ls will show 
# raw SLC intf topo
#
alias rm 'rm -f'
unset noclobber
#
if ($#argv < 1) then
 echo " "
 echo "Usage: cleanup.csh directory "
 echo " "
 echo " directory could be: raw, SLC, topo, intf, or all "
 echo " "
 echo "Example: cleanup.csh all"
 echo " "
 exit 1
endif
#
#
if( $1 == all) then
  echo ""
  echo "clean up all"
  rm -r SLC
  rm -r intf
  rm raw/*.PRM*
  rm raw/*.raw 
  rm raw/*.LED
  rm raw/*.SLC
  cd topo
  ls | grep -v dem.grd | xargs rm
  cd ..
  echo ""
endif
if( $1 == raw) then
  echo ""
  echo "clean up raw/ folder"
  rm raw/*.PRM*
  rm raw/*.raw
  rm raw/*.LED
  echo ""
endif
if( $1 == SLC) then
  echo ""
  echo "clean up SLC/ folder"
  rm -r SLC/*
  echo ""
endif
if( $1 == intf) then
  echo ""
  echo "clean up intf/ folder"
  rm -r intf/*
  echo ""
endif
if ( $1 == topo ) then 
  echo ""
  echo "clean up topo/ folder"
  cd topo
  ls | grep -v dem.grd | xargs rm 
  cd ..
  echo ""
endif
