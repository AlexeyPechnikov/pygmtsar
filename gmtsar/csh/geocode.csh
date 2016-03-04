#!/bin/csh -f
#       $Id$
#
#  D. Sandwell FEB 10 2010
#  Kurt Feigl 20150811 add annotation to grd files
#
alias rm 'rm -f'
unset noclobber
#
  if ($#argv < 1) then
errormessage:
    echo ""
    echo "Usage: geocode.csh correlation_threshold"
    echo ""
    echo " phase is masked when correlation is less than correlation_threshold"
    echo ""
    echo "Example: geocode.csh .12"
    echo ""
    exit 1
  endif
#
#   first mask the phase and phase gradient using the correlation
#
gmt grdmath corr.grd $1 GE 0 NAN mask.grd MUL = mask2.grd -V
gmt grdmath phase.grd mask2.grd MUL = phase_mask.grd
if (-e xphase.grd) then
  gmt grdmath xphase.grd mask2.grd MUL = xphase_mask.grd
  gmt grdmath yphase.grd mask2.grd MUL = yphase_mask.grd
endif
if (-e unwrap.grd) then 
  gmt grdcut mask2.grd `gmt grdinfo unwrap.grd -I-` -Gmask3.grd
  gmt grdmath unwrap.grd mask3.grd MUL = unwrap_mask.grd
endif
if (-e phasefilt.grd) then 
  gmt grdmath phasefilt.grd mask2.grd MUL = phasefilt_mask.grd
endif
#
#   look at the masked phase
#
set boundR = `gmt grdinfo display_amp.grd -C | awk '{print ($3-$2)/4}'`
set boundA = `gmt grdinfo display_amp.grd -C | awk '{print ($5-$4)/4}'`
gmt grdimage phase_mask.grd -JX6.5i -Cphase.cpt -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3i -Y3i -P -K > phase_mask.ps
gmt psscale -D3.3/-1.5/5/0.2h -Cphase.cpt -B1.57:"phase, rad": -O >> phase_mask.ps
if (-e xphase_mask.grd) then
  gmt grdimage xphase_mask.grd -JX8i -Cphase_grad.cpt -X.2i -Y.5i -P > xphase_mask.ps
  gmt grdimage yphase_mask.grd -JX8i -Cphase_grad.cpt -X.2i -Y.5i -P > yphase_mask.ps
endif
if (-e unwrap_mask.grd) then 
  gmt grdimage unwrap_mask.grd -JX6.5i -B"$boundR":Range:/"$boundA":Azimuth:WSen -Cunwrap.cpt -X1.3i -Y3i -P -K > unwrap_mask.ps
  set std = `gmt grdinfo -C -L2 unwrap_mask.grd | awk '{printf("%5.1f", $13)}'`
  gmt psscale -D3.3/-1.5/5/0.2h -Cunwrap.cpt -B"$std":"unwrapped phase, rad": -O -E >> unwrap_mask.ps
endif
if (-e phasefilt_mask.grd) then 
  gmt grdimage phasefilt_mask.grd -JX6.5i -B"$boundR":Range:/"$boundA":Azimuth:WSen -Cphase.cpt -X1.3i -Y3i -P -K > phasefilt_mask.ps
  gmt psscale -D3.3/-1.5/5/0.2h -Cphase.cpt -B1.57:"phase, rad": -O >> phasefilt_mask.ps
endif
# line-of-sight displacement
if (-e unwrap_mask.grd) then
  set wavel = `grep wavelength *.PRM | awk '{print($3)}' | head -1 `
  gmt grdmath unwrap_mask.grd $wavel MUL -79.58 MUL = los.grd
  gmt grdgradient los.grd -Nt.9 -A0. -Glos_grad.grd
  set tmp = `gmt grdinfo -C -L2 los.grd`
  set limitU = `echo $tmp | awk '{printf("%5.1f", $12+$13*2)}'`
  set limitL = `echo $tmp | awk '{printf("%5.1f", $12-$13*2)}'`
  set std = `echo $tmp | awk '{printf("%5.1f", $13)}'`
  gmt makecpt -Cpolar -Z -T"$limitL"/"$limitU"/1 -D > los.cpt
  gmt grdimage los.grd -Ilos_grad.grd -Clos.cpt -B"$boundR":Range:/"$boundA":Azimuth:WSen -JX6.5i -X1.3i -Y3i -P -K > los.ps
  gmt psscale -D3.3/-1.5/4/0.2h -Clos.cpt -B"$std":"LOS displacement, mm":/:"range decrease": -O -E >> los.ps 
endif

#
#  now reproject the phase to lon/lat space
#
echo "geocode.csh"
echo "project correlation, phase, unwrapped and amplitude back to lon lat coordinates"
set maker = $0:t
set today = `date`
set remarked = `echo by $USER on $today with $maker`
echo remarked is $remarked

 proj_ra2ll.csh trans.dat corr.grd        corr_ll.grd           ; gmt grdedit -D//"dimensionless"/1///"$PWD:t geocoded correlation"/"$remarked"      corr_ll.grd
#proj_ra2ll.csh trans.dat phase.grd       phase_ll.grd          ; gmt grdedit -D//"radians"/1///"$PWD:t wrapped phase"/"$remarked"                   phase_ll.grd
 proj_ra2ll.csh trans.dat phasefilt.grd   phasefilt_ll.grd      ; gmt grdedit -D//"radians"/1///"$PWD:t wrapped phase after filtering"/"$remarked"   phasefilt_ll.grd
proj_ra2ll.csh trans.dat phase_mask.grd  phase_mask_ll.grd     ; gmt grdedit -D//"radians"/1///"$PWD:t wrapped phase after masking"/"$remarked"     phase_mask_ll.grd
 proj_ra2ll.csh trans.dat display_amp.grd display_amp_ll.grd    ; gmt grdedit -D//"dimensionless"/1///"PWD:t amplitude"/"$remarked"                  display_amp_ll.grd
if (-e xphase_mask.grd) then
  proj_ra2ll.csh trans.dat xphase_mask.grd xphase_mask_ll.grd  ; gmt grdedit -D//"radians"/1///PWD:t xphase"/"$remarked"                            xphase_mask_ll.grd
  proj_ra2ll.csh trans.dat yphase_mask.grd yphase_mask_ll.grd  ; gmt grdedit -D//"radians"/1///PWD:t yphase"/"$remarked"                            yphase_mask_ll.grd
endif
if (-e unwrap_mask.grd) then
  proj_ra2ll.csh trans.dat unwrap_mask.grd unwrap_mask_ll.grd  ; gmt grdedit -D//"radians"/1///"PWD:t unwrapped, masked phase"/"$remarked"               unwrap_mask_ll.grd
endif
if (-e unwrap.grd) then
  proj_ra2ll.csh trans.dat unwrap.grd unwrap_ll.grd  ; gmt grdedit -D//"radians"/1///"PWD:t unwrapped phase"/"$remarked"               unwrap_ll.grd
endif
if (-e phasefilt_mask.grd) then
  proj_ra2ll.csh trans.dat phasefilt_mask.grd phasefilt_mask_ll.grd ; gmt grdedit -D//"phase in radians"/1///"PWD:t wrapped phase masked filtered"/"$remarked"   phasefilt_mask_ll.grd
endif


#
#   now image for google earth
#
echo "geocode.csh"
echo "make the KML files for Google Earth"
grd2kml.csh display_amp_ll display_amp.cpt
grd2kml.csh corr_ll corr.cpt
grd2kml.csh phase_mask_ll phase.cpt
grd2kml.csh phasefilt_mask_ll phase.cpt
#ln -s phasefilt_mask_ll.grd phase_mask_ll_bw.grd
#grd2kml.csh phase_mask_ll_bw phase_bw.cpt
#rm phase_mask_ll_bw.grd
if (-e xphase_mask_ll.grd) then
  grd2kml.csh xphase_mask_ll phase_grad.cpt
  grd2kml.csh yphase_mask_ll phase_grad.cpt
endif
if (-e unwrap_mask_ll.grd) then
  grd2kml.csh unwrap_mask_ll unwrap.cpt
endif
if (-e phasefilt_mask_ll.grd) then
  grd2kml.csh phasefilt_mask_ll phase.cpt
endif
if (-e unwrap_mask_ll.grd) then
  # constant is negative to make LOS = -1 * range change
  # constant is (1000 mm) / (4 * pi)
   gmt grdmath unwrap_mask_ll.grd $wavel MUL -79.58 MUL = los_ll.grd 

   gmt grdedit -D//"mm"/1///"$PWD:t LOS displacement"/"equals negative range" los_ll.grd 

  grd2kml.csh los_ll los.cpt
endif
