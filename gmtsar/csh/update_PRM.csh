#!/bin/csh -f
#       $Id$
#
#  D. T. Sandwell AUG 30 2007
#
#  Script to replace a value in a PRM-file
#
#
alias rm 'rm -f'
unset noclobber
#
 if ($#argv < 2) then
  echo " "
  echo "Usage: update_PRM file.PRM param value "
  echo " "
  echo "Example: update_PRM IMG-HH-ALPSRP049040660-H1.0__A.PRM rshift 10"
  exit 1
 endif
#
#  copy the file to junk.PRM
#
 cp $1 junk.PRM
#
#  remove the line with the matching parameter
#
 if ($2 != 'clock_start' && $2 != 'clock_stop') then
   grep -v $2 < junk.PRM > junk2.PRM
 else
   if ($2 == 'clock_start') then
     set tmptime = `grep SC_clock_start < junk.PRM | awk '{print $3}'`
     grep -v $2 < junk.PRM > junk2.PRM
     echo 'SC_clock_start' ' = ' $tmptime >> junk2.PRM
   else
     set tmptime = `grep SC_clock_stop < junk.PRM | awk '{print $3}'`
     grep -v $2 < junk.PRM > junk2.PRM
     echo 'SC_clock_stop' ' = ' $tmptime >> junk2.PRM
   endif
 endif
#
#  add a new line
#
 echo $2 ' = ' $3 >> junk2.PRM
#
 mv junk2.PRM $1
 rm junk.PRM
