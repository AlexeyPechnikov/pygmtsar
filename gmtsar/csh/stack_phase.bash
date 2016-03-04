#!/bin/bash -f

#  compute the mean LOS velocity and its standard deviations from a stack of unwrapped phase

# Creator: Xiaopeng Tong, David Sandwell
# Date:    July 6th 2013

# Reference: GMTSAR document Appendix C

if [ $# -ne 3 ]; then
  echo ""
  echo "Usage: stack_phase.bash stack.list mean.grd std.grd"
  echo ""
  echo "  compute the mean LOS velocity (mm/yr) and its standard deviations"
  echo "  from a stack of unwrapped phase" 
  echo ""
  echo "  Note the unwrapped phase can be pre-processed (detrended/filtered, etc)"
  echo ""
  echo "    stack.list  --     phase grd files to be stacked "
  echo "                       with PRM files corresponding to the interferogram"
  echo ""
  echo "    example of the stack.list:"
  echo "     unwrap_03562_05575.grd IMG-HH-ALPSRP035620660-H1.0__A.PRM IMG-HH-ALPSRP055750660-H1.0__A.PRM"
  echo "     unwrap_03562_12956.grd IMG-HH-ALPSRP035620660-H1.0__A.PRM IMG-HH-ALPSRP129560660-H1.0__A.PRM" 
  echo "     unwrap_05575_12956.grd IMG-HH-ALPSRP055750660-H1.0__A.PRM IMG-HH-ALPSRP129560660-H1.0__A.PRM"
  echo ""
  echo "  Make sure the required grd files and PRM files exist"
  echo ""
  echo "    mean.grd    --     output file: mean LOS velocity"
  echo "    std.grd     --     output file: standard deviation from the mean"
  echo ""
  exit 1
fi

if [ ! -e $1 ] 
then
  echo ""
  echo "no input file found: $1"
  echo ""
  exit 1
fi

list=$1
outmean=$2
outstd=$3

# compute the mean LOS velocity from the unwrapped phase
echo 
echo " compute the mean LOS velocity .. "
echo 

let "ac=0"
exec<$list
declare -a dayarray=()

while read grd ref rep 
do
  # echo $grd $ref $rep

  if [ ! -e $grd ];
  then
    echo "file does not exist: $grd"
    exit 1
  fi

  year_ref=`grep SC_clock_start $ref | awk '{print $3}' | cut -c1-4`
  day_ref=`grep SC_clock_start $ref | awk '{print $3}' | cut -c5-7`
  year_rep=`grep SC_clock_start $rep | awk '{print $3}' | cut -c1-4`
  day_rep=`grep SC_clock_start $rep | awk '{print $3}' | cut -c5-7`

  numdays=`echo "($year_rep-$year_ref)*365+$day_rep-$day_ref" | bc`

  dayarray=(${dayarray[@]} $numdays)  

  if [ "$numdays" -lt 0 ]; then
    echo "time span should be >= 0"
    exit 1
  else
    echo "time span (days) of $grd:" $numdays
  fi

  if [ "$ac" -eq 0 ]; then
    gmt grdmath $grd = sum.grd 
    wavel=`grep wavelength $ref | awk '{print($3)}'`
  else
    gmt grdmath $grd sum.grd ADD = tmp.grd 
    mv tmp.grd sum.grd  
  fi

  let "ac += $numdays"
  echo "accumulative days in the stack:" $ac
done
gmt grdmath sum.grd $ac DIV 365 MUL $wavel MUL -79.58 MUL = $outmean

# compute the standard deviation from the mean LOS velocity 
echo 
echo " compute the standard deviations .. "
echo 

let "num=0"
exec<$list
while read grd ref rep
do
  if [ "$num" -eq 0 ]; then
    gmt grdmath $grd ${dayarray[$num]} DIV 365 MUL $wavel MUL -79.58 MUL $outmean SUB SQR = sum2.grd 
  else 
    gmt grdmath $grd ${dayarray[$num]} DIV 365 MUL $wavel MUL -79.58 MUL $outmean SUB SQR sum2.grd ADD = sum2tmp.grd
    mv sum2tmp.grd sum2.grd
  fi
  let "num += 1"
done 
gmt grdmath sum2.grd $num DIV SQRT = $outstd

# clean up
rm sum.grd sum2.grd 
