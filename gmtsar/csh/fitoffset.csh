#!/bin/csh -f
#
#  D. Sandwell DEC 10 2007
#
#  compute 2-6 alignment parameters from xcorr.dat
#
alias rm 'rm -f'
unset noclobber
#
# check for number of arguments
#
if ($#argv == 0) then
  	echo "  "
  	echo "Usage: fitoffset.csh npar xcorr.dat [SNR]"
  	echo "  "
        echo "        npar        - number of parameters to fit  "
  	echo "        xcorr.dat   - files of range and azimuth offset estimates "
  	echo "        SNR         - optional SNR cutoff (default 20)"
  	echo "  "
  	echo "Example: fitoffset.csh 6 freq_xcorr.dat "
  	echo "  "
  	exit 1
endif
#
 if ($#argv == 3) then
  	set SNR = $3
 else
  	set SNR = 20.
 	endif
endif
#
#  first extract the range and azimuth data
#
 awk '{ if ($5 > '$SNR') printf("%f %f %f \n",$1,$3,$2); }' < $2 > r.xyz
 awk '{ if ($5 > '$SNR') printf("%f %f %f \n",$1,$3,$4); }' < $2 > a.xyz
#
#  make sure there are enough points remaining, otherwise exit
#
 set NPTS0 = ` wc -l $2 | awk '{print $1}'`
 set NPTS = ` wc -l r.xyz | awk '{print $1}'`
 if($NPTS < 8) then
  echo "  "
  echo " FAILED - not enough points to estimate parameters"
  echo " try lower SNR "
  echo " NPTS0, NPTS  " $NPTS0 $NPTS
  echo "  "
  exit 1
 endif
#
# compute requested number of parameters 2-6
#
 if($1 == 2) then
#
#   use gmt trend2d
#
 ( gmt trend2d r.xyz -Fxyz -N1r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > rm.coef
 ( gmt trend2d a.xyz -Fxyz -N1r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > am.coef
#
#  get the data range and paste to the coeffifients
#
 gmt gmtinfo r.xyz -C | awk '{print $1, $2, $3, $4 }' > range.coef
 paste rm.coef range.coef  > rm.rng
 paste am.coef range.coef  > am.rng
 rm r.xyz a.xyz am.coef rm.coef range.coef
#
#   now convert to range coefficients
#
 awk '{print ( $1 ) }' < rm.rng > rshift
 awk '{if($1 >= 0) {printf("%s %g \n","rshift = ",int($1)); printf("%s %g \n","sub_int_r = ",($1 %1))}}' < rshift
 awk '{if($1 < 0)  {printf("%s %g \n","rshift = ",int($1)-1); printf("%s %g \n","sub_int_r = ",($1 %1)+1)}}' < rshift
 awk '{printf ("%s %g \n","stretch_r = ",0.0)}' < rm.rng
 awk '{printf ("%s %g \n","a_stretch_r = ",0.0)}' < rm.rng
#
#  now convert to azimuth coefficients
#
 awk '{print ( $1 ) }' < am.rng > ashift
 awk '{if($1 >= 0) {printf("%s %g \n","ashift = ",int($1)); printf("%s %g \n","sub_int_a = ",($1 %1))}}' < ashift
 awk '{if($1 < 0)  {printf("%s %g \n","ashift = ",int($1)-1); printf("%s %g \n","sub_int_a = ",($1 %1)+1)}}' < ashift
 awk '{printf ("%s %g \n","stretch_a = ",0.0)}' < am.rng
 awk '{printf ("%s %g \n","a_stretch_a = ",0.0)}' < am.rng
#
 else if ($1 == 3) then
#
#   use gmt trend2d
#
 ( gmt trend2d r.xyz -Fxyz -N2r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > rm.coef
 ( gmt trend2d a.xyz -Fxyz -N1r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > am.coef
#
#  get the data range and paste to the coeffifients
#
 gmt gmtinfo r.xyz -C | awk '{print $1, $2, $3, $4 }' > range.coef
 paste rm.coef range.coef  > rm.rng
 paste am.coef range.coef  > am.rng
 rm r.xyz a.xyz am.coef rm.coef range.coef
#
#   now convert to range coefficients
#
 awk '{print ( $1 - $2*($4+$3)/($4-$3) ) }' < rm.rng > rshift
 awk '{if($1 >= 0) {printf("%s %g \n","rshift = ",int($1)); printf("%s %g \n","sub_int_r = ",($1 %1))}}' < rshift
 awk '{if($1 < 0)  {printf("%s %g \n","rshift = ",int($1)-1); printf("%s %g \n","sub_int_r = ",($1 %1)+1)}}' < rshift
 awk '{printf ("%s %g \n","stretch_r = ",$2*2./($4-$3))}' < rm.rng
 awk '{printf ("%s %g \n","a_stretch_r = ",0.0)}' < rm.rng
#
#  now convert to azimuth coefficients
#
 awk '{print ( $1 ) }' < am.rng > ashift
 awk '{if($1 >= 0) {printf("%s %g \n","ashift = ",int($1)); printf("%s %g \n","sub_int_a = ",($1 %1))}}' < ashift
 awk '{if($1 < 0)  {printf("%s %g \n","ashift = ",int($1)-1); printf("%s %g \n","sub_int_a = ",($1 %1)+1)}}' < ashift
 awk '{printf ("%s %g \n","stretch_a = ",0.0)}' < am.rng
 awk '{printf ("%s %g \n","a_stretch_a = ",0.0)}' < am.rng
#
 else if ($1 == 4) then
#
#   use gmt trend2d
#
 ( gmt trend2d r.xyz -Fxyz -N2r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > rm.coef
 ( gmt trend2d a.xyz -Fxyz -N2r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > am.coef
#
#  get the data range and paste to the coeffifients
#
 gmt gmtinfo r.xyz -C | awk '{print $1, $2, $3, $4 }' > range.coef
 paste rm.coef range.coef  > rm.rng
 paste am.coef range.coef  > am.rng
 rm r.xyz a.xyz am.coef rm.coef range.coef
#
#   now convert to range coefficients
#
 awk '{print ( $1 - $2*($4+$3)/($4-$3) ) }' < rm.rng > rshift
 awk '{if($1 >= 0) {printf("%s %g \n","rshift = ",int($1)); printf("%s %g \n","sub_int_r = ",($1 %1))}}' < rshift
 awk '{if($1 < 0)  {printf("%s %g \n","rshift = ",int($1)-1); printf("%s %g \n","sub_int_r = ",($1 %1)+1)}}' < rshift
 awk '{printf ("%s %g \n","stretch_r = ",$2*2./($4-$3))}' < rm.rng
 awk '{printf ("%s %g \n","a_stretch_r = ",0.0)}' < rm.rng
#
#  now convert to azimuth coefficients
#
 awk '{print ( $1 - $2*($4+$3)/($4-$3) ) }' < am.rng > ashift
 awk '{if($1 >= 0) {printf("%s %g \n","ashift = ",int($1)); printf("%s %g \n","sub_int_a = ",($1 %1))}}' < ashift
 awk '{if($1 < 0)  {printf("%s %g \n","ashift = ",int($1)-1); printf("%s %g \n","sub_int_a = ",($1 %1)+1)}}' < ashift
 awk '{printf ("%s %g \n","stretch_a = ",$2*2./($4-$3))}' < am.rng
 awk '{printf ("%s %g \n","a_stretch_a = ",0.0)}' < am.rng
#
 else if ($1 == 5) then
#
#   use gmt trend2d
#
 ( gmt trend2d r.xyz -Fxyz -N2r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > rm.coef
 ( gmt trend2d a.xyz -Fxyz -N3r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > am.coef
#
#  get the data range and paste to the coeffifients
#
 gmt gmtinfo r.xyz -C | awk '{print $1, $2, $3, $4 }' > range.coef
 paste rm.coef range.coef  > rm.rng
 paste am.coef range.coef  > am.rng
 rm r.xyz a.xyz am.coef rm.coef range.coef
#
#   now convert to range coefficients
#
 awk '{print ( $1 - $2*($4+$3)/($4-$3) ) }' < rm.rng > rshift
 awk '{if($1 >= 0) {printf("%s %g \n","rshift = ",int($1)); printf("%s %g \n","sub_int_r = ",($1 %1))}}' < rshift
 awk '{if($1 < 0)  {printf("%s %g \n","rshift = ",int($1)-1); printf("%s %g \n","sub_int_r = ",($1 %1)+1)}}' < rshift
 awk '{printf ("%s %g \n","stretch_r = ",$2*2./($4-$3))}' < rm.rng
 awk '{printf ("%s %g \n","a_stretch_r = ",0.0)}' < rm.rng
#
#  now convert to azimuth coefficients
#
 awk '{print ( $1 - $2*($5+$4)/($5-$4) -$3*($7+$6)/($7-$6) ) }' < am.rng > ashift
 awk '{if($1 >= 0) {printf("%s %g \n","ashift = ",int($1)); printf("%s %g \n","sub_int_a = ",($1 %1))}}' < ashift
 awk '{if($1 < 0)  {printf("%s %g \n","ashift = ",int($1)-1); printf("%s %g \n","sub_int_a = ",($1 %1)+1)}}' < ashift
 awk '{printf ("%s %g \n","stretch_a = ",$2*2./($5-$4))}' < am.rng
 awk '{printf ("%s %g \n","a_stretch_a = ",$3*2./($7-$6))}' < am.rng
#
 else if ($1 == 6) then
#
#   use gmt trend2d
#
 ( gmt trend2d r.xyz -Fxyz -N3r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > rm.coef
 ( gmt trend2d a.xyz -Fxyz -N3r -V >  /dev/null ) |& grep oefficients | awk -F":" '{print $NF'} > am.coef
#
#  get the data range and paste to the coeffifients
#
 gmt gmtinfo r.xyz -C | awk '{print $1, $2, $3, $4 }' > range.coef
 paste rm.coef range.coef  > rm.rng
 paste am.coef range.coef  > am.rng
 rm r.xyz a.xyz am.coef rm.coef range.coef
#
#   now convert to range coefficients
#
 awk '{print ( $1 - $2*($5+$4)/($5-$4) -$3*($7+$6)/($7-$6) ) }' < rm.rng > rshift
 awk '{if($1 >= 0) {printf("%s %g \n","rshift = ",int($1)); printf("%s %g \n","sub_int_r = ",($1 %1))}}' < rshift
 awk '{if($1 < 0)  {printf("%s %g \n","rshift = ",int($1)-1); printf("%s %g \n","sub_int_r = ",($1 %1)+1)}}' < rshift
 awk '{printf ("%s %g \n","stretch_r = ",$2*2./($5-$4))}' < rm.rng
 awk '{printf ("%s %g \n","a_stretch_r = ",$3*2./($7-$6))}' < rm.rng
#
#  now convert to azimuth coefficients
#
 awk '{print ( $1 - $2*($5+$4)/($5-$4) -$3*($7+$6)/($7-$6) ) }' < am.rng > ashift
 awk '{if($1 >= 0) {printf("%s %g \n","ashift = ",int($1)); printf("%s %g \n","sub_int_a = ",($1 %1))}}' < ashift
 awk '{if($1 < 0)  {printf("%s %g \n","ashift = ",int($1)-1); printf("%s %g \n","sub_int_a = ",($1 %1)+1)}}' < ashift
 awk '{printf ("%s %g \n","stretch_a = ",$2*2./($5-$4))}' < am.rng
 awk '{printf ("%s %g \n","a_stretch_a = ",$3*2./($7-$6))}' < am.rng
#
 else
  echo "number of parameters must be either 2,3,4,5 or 6 "
 endif
#
#  cleanup
#
 rm rshift ashift rm.rng am.rng
#
#  OK we are done
#
#  echo " NPTS0, NPTS  " $NPTS0 $NPTS
