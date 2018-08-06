#!/bin/csh -f
#       $Id$
#
#  D. Sandwell FEB 7 2010
#
unset noclobber
#
#  script to focus SAR data and deal with negative yshift
#
# check for number of arguments
#
 if ($#argv < 1) then
  echo ""
  echo "Usage: sarp.csh file.PRM "
  echo ""
  echo "Example: sarp.csh IMG-HH-ALPSRP049040660-H1.0__A.PRM "
  echo ""
  exit 1
 endif
#
 set  slcfile = ` echo $1 | awk '{ print(substr($1,1,length($1)-3)"SLC")}'`
#
 esarp  $1 $slcfile > tmp_sarp
#
#   update the PRM file
#
grep -v SLC_file $1 | grep -v dtype > tmp.PRM
mv tmp.PRM $1
echo "SLC_file  = "$slcfile >> $1
if ( $3 == R4) then
	echo "dtype = c" >> $1
else
	echo "dtype = a" >> $1
endif

# 
#   put in dfact into PRM file
#
set dfact = `grep I2SCALE tmp_sarp | awk '{print $2}'`
update_PRM $1 SLC_scale $dfact
rm tmp_sarp
