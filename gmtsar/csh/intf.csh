#!/bin/csh -f
#       $Id$
#
#  Xiaopeng Tong FEB 4 2010
# 
# Make the interferogram.
#
#  Matt Wei May 4 2010, ENVISAT
#
# add in TSX, Jan 10 2014
alias rm 'rm -f'
gmt set IO_NC4_CHUNK_SIZE classic
#
#
  if ($#argv < 2) then
errormessage:
    echo ""
    echo "Usage: intf.csh ref.PRM rep.PRM [-topo topogrd] [-model modelgrd]"
    echo ""
    echo " The dimensions of topogrd and modelgrd should be compatible with SLC file."
    echo ""
    echo "Example: intf.csh IMG-HH-ALPSRP055750660-H1.0__A.PRM IMG-HH-ALPSRP049040660-H1.0__A.PRM -topo topo_ra.grd"
    echo ""
    exit 1
  endif
#
# which satellite
#
  set SC = `grep SC_identity $1 | awk '{print $3}'`
  if ($SC == 1 || $SC == 2 || $SC == 4 || $SC == 5 || $SC == 6) then
    cp $2 $2"0"
    cp $1 $1"0"
    SAT_baseline $1 $2 | tail -n6 >> $2
    SAT_baseline $1 $1 | grep height >> $1
  else if ($SC > 6) then
    cp $2 $2"0"
    cp $1 $1"0"
    SAT_baseline $1 $2 | tail -n6 >> $2
    SAT_baseline $1 $1 | grep height >> $1
  else
    echo "Incorrect satellite id in prm file"
    exit 0
  endif
#
# form the interferogram optionally using topo_ra and modelphase
# 
  if ($#argv == 2 || $#argv == 4 || $#argv == 6) then 
    echo "intf.csh"
    echo "running phasediff..."
    phasediff $argv 
  else
    goto errormessage
  endif
  mv $1"0" $1
  mv $2"0" $2
