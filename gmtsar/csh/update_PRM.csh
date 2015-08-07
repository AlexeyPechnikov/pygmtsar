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
 grep -v $2 < junk.PRM > junk2.PRM
#
#  add a new line
#
 echo $2 ' = ' $3 >> junk2.PRM
#
 mv junk2.PRM $1
 rm junk.PRM
