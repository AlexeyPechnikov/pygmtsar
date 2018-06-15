#!/bin/csh -f
#       $Id$
#
#  Xiaopeng Tong and David Sandwell 
#  FEB 4 2010
#  Matt Wei May 4 2010, ENVISAT
#  DTS - May 26, 2010, added phase gadient
#  EF, DTS, XT - Jan 10 2014, TSX
#
# Convolve the real.grd and imag.grd with gaussian filters. 
# Form amplitude, phase, phase gradient, and correlation images. 
#
#
  alias rm 'rm -f'
  gmt set IO_NC4_CHUNK_SIZE classic
#
#
# set grdimage options
#
  set scale = "-JX6.5i"
  set thresh = "5.e-21"
  gmt set COLOR_MODEL = hsv
  gmt set PROJ_LENGTH_UNIT = inch

  if ($#argv != 4 && $#argv != 6) then
errormessage:
    echo ""
    echo "Usage: filter.csh master.PRM slave.PRM filter decimation [rng_dec azi_dec]"
    echo ""
    echo " Apply gaussian filter to amplitude and phase images."
    echo " "
    echo " filter -  wavelength of the filter in meters (0.5 gain)"
    echo " decimation - (1) better resolution, (2) smaller files"
    echo " "
    echo "Example: filter.csh IMG-HH-ALPSRP055750660-H1.0__A.PRM IMG-HH-ALPSRP049040660-H1.0__A.PRM 300  2"
    echo ""
    exit 1
  endif
  echo "filter.csh"
#
# define filter and decimation variables
#
  set sharedir = `gmtsar_sharedir.csh`
  set filter3 = $sharedir/filters/fill.3x3
  set filter4 = $sharedir/filters/xdir
  set filter5 = $sharedir/filters/ydir
  set dec  = $4
  set az_lks = 4 
  set PRF = `grep PRF *.PRM | awk 'NR == 1 {printf("%d", $3)}'`
  if( $PRF < 1000 ) then
     set az_lks = 1
  endif
#
# look for range sampling rate
#
  set rng_samp_rate = `grep rng_samp_rate $1 | awk 'NR == 1 {printf("%d", $3)}'`
#
# set the range spacing in units of image range pixel size
#
  if ($?rng_samp_rate) then
    if ($rng_samp_rate > 110000000) then 
      set dec_rng = 4
      set filter1 = $sharedir/filters/gauss15x5
    else if ($rng_samp_rate < 110000000 && $rng_samp_rate > 20000000) then
      set dec_rng = 2
      set filter1 = $sharedir/filters/gauss15x5
#
# special for TOPS mode
#
      if($az_lks == 1) then
        set filter1 = $sharedir/filters/gauss5x5
      endif
    else  
      set dec_rng = 1
      set filter1 = $sharedir/filters/gauss15x3
    endif
  else
    echo "Undefined rng_samp_rate in the master PRM file"
    exit 1
  endif
#
#  make the custom filter2 and set the decimation
#
  make_gaussian_filter $1 $dec_rng $az_lks $3 > ijdec
  set filter2 = gauss_$3
  set idec = `cat ijdec | awk -v dc="$dec" '{ print dc*$1 }'`
  set jdec = `cat ijdec | awk -v dc="$dec" '{ print dc*$2 }'`
  if($#argv == 6) then
    set idec = `echo $6 $az_lks | awk '{printf("%d",$1/$2)}'`
    set jdec = `echo $5 $dec_rng | awk '{printf("%d",$1/$2)}'`
    echo "setting range_dec = $5, azimuth_dec = $6"
  endif
  echo "$filter2 $idec $jdec ($az_lks $dec_rng)" 
#
# filter the two amplitude images
#
  echo "making amplitudes..."
  conv $az_lks $dec_rng $filter1 $1 amp1_tmp.grd=bf
  conv $idec $jdec $filter2 amp1_tmp.grd=bf amp1.grd
  rm amp1_tmp.grd
  conv $az_lks $dec_rng $filter1 $2 amp2_tmp.grd=bf
  conv $idec $jdec $filter2 amp2_tmp.grd=bf amp2.grd
  rm amp2_tmp.grd
#
# filter the real and imaginary parts of the interferogram
# also compute gradients
#
  echo "filtering interferogram..."
  conv $az_lks $dec_rng $filter1 real.grd=bf real_tmp.grd=bf
  conv $idec $jdec $filter2 real_tmp.grd=bf realfilt.grd
#  conv $dec $dec $filter4 real_tmp.grd xreal.grd
#  conv $dec $dec $filter5 real_tmp.grd yreal.grd
  rm real_tmp.grd 
  rm real.grd
  conv $az_lks $dec_rng $filter1 imag.grd=bf imag_tmp.grd=bf
  conv $idec $jdec $filter2 imag_tmp.grd=bf imagfilt.grd
#  conv $dec $dec $filter4 imag_tmp.grd ximag.grd
#  conv $dec $dec $filter5 imag_tmp.grd yimag.grd
  rm imag_tmp.grd 
  rm imag.grd
#
# form amplitude image
#
  echo "making amplitude..."
  gmt grdmath realfilt.grd imagfilt.grd HYPOT  = amp.grd 
  gmt grdmath amp.grd 0.5 POW FLIPUD = display_amp.grd 
  set AMAX = `gmt grdinfo -L2 display_amp.grd | grep stdev | awk '{ print 3*$5 }'`
  gmt grd2cpt display_amp.grd -Z -D -L0/$AMAX -Cgray > display_amp.cpt
  echo "N  255   255   254" >> display_amp.cpt
  gmt grdimage display_amp.grd -Cdisplay_amp.cpt $scale -Bxaf+lRange -Byaf+lAzimuth -BWSen -X1.3i -Y3i -P -K > display_amp.ps
  gmt psscale -Rdisplay_amp.grd -J -DJTC+w5i/0.2i+h+ef -Cdisplay_amp.cpt -Bx0+l"Amplitude (histogram equalized)" -O >> display_amp.ps
  gmt psconvert -Tf -P -Z display_amp.ps
  echo "Amplitude map: display_amp.pdf"
#
# form the correlation
#
  echo "making correlation..."
  gmt grdmath amp1.grd amp2.grd MUL = tmp.grd
  gmt grdmath tmp.grd $thresh GE 0 NAN = mask.grd
  gmt grdmath amp.grd tmp.grd SQRT DIV mask.grd MUL FLIPUD = tmp2.grd=bf
  conv 1 1 $filter3 tmp2.grd=bf corr.grd
  gmt makecpt -T0./.8/0.1 -Cgray -Z -N > corr.cpt
  echo "N  255   255   254" >> corr.cpt
  gmt grdimage corr.grd $scale -Ccorr.cpt -Bxaf+lRange -Byaf+lAzimuth -BWSen -X1.3i -Y3i -P -K > corr.ps
  gmt psscale -Rcorr.grd -J -DJTC+w5i/0.2i+h+ef -Ccorr.cpt -Baf+lCorrelation -O >> corr.ps
  gmt psconvert -Tf -P -Z corr.ps
  echo "Correlation map: corr.pdf"
#
# form the phase 
#
  echo "making phase..."
  gmt grdmath imagfilt.grd realfilt.grd ATAN2 mask.grd MUL FLIPUD = phase.grd
  gmt makecpt -Crainbow -T-3.15/3.15/0.1 -Z -N > phase.cpt
  gmt grdimage phase.grd $scale -Bxaf+lRange -Byaf+lAzimuth -BWSen -Cphase.cpt -X1.3i -Y3i -P -K > phase.ps
  gmt psscale -Rphase.grd -J -DJTC+w5i/0.2i+h -Cphase.cpt -B1.57+l"Phase" -By+lrad -O >> phase.ps
  gmt psconvert -Tf -P -Z phase.ps
  echo "Phase map: phase.pdf"
#
# compute the solid earth tide
# uncomment lines with ##
#
##ln -s ../../topo/dem.grd .
##tide_correction.csh $1 $2 dem.grd
##mv tide.grd tmp.grd
##gmt grdsample tmp.grd -Rimagfilt.grd -Gtide.grd
#
#  make the Werner/Goldstein filtered phase
#
  echo "filtering phase..."
  phasefilt -imag imagfilt.grd -real realfilt.grd -amp1 amp1.grd -amp2 amp2.grd -psize 32 
  gmt grdedit filtphase.grd `gmt grdinfo mask.grd -I- --FORMAT_FLOAT_OUT=%.12lg` 
  gmt grdmath filtphase.grd mask.grd MUL FLIPUD = phasefilt.grd
##cp phasefilt.grd phasefilt_old.grd
##gmt grdmath phasefilt.grd tide.grd SUB PI ADD 2 PI MUL MOD PI SUB = phasefilt.grd
  rm filtphase.grd
  gmt grdimage phasefilt.grd $scale -Bxaf+lRange -Byaf+lAzimuth -BWSen -Cphase.cpt -X1.3i -Y3i -P -K > phasefilt.ps
  gmt psscale -Rphasefilt.grd -J -DJTC+w5i/0.2i+h -Cphase.cpt -Bxa1.57+l"Phase" -By+lrad -O >> phasefilt.ps
  gmt psconvert -Tf -P -Z phasefilt.ps
  echo "Filtered phase map: phasefilt.pdf"
# 
#  form the phase gradients
#
#  echo "making phase gradient..."
#  gmt grdmath amp.grd 2. POW = amp_pow.grd
#  gmt grdmath realfilt.grd ximag.grd MUL imagfilt.grd xreal.grd MUL SUB amp_pow.grd DIV mask.grd MUL FLIPUD = xphase.grd
#  gmt grdmath realfilt.grd yimag.grd MUL imagfilt.grd yreal.grd MUL SUB amp_pow.grd DIV mask.grd MUL FLIPUD = yphase.grd 
#  gmt makecpt -Cgray -T-0.7/0.7/0.1 -Z -N > phase_grad.cpt
#  echo "N  255   255   254" >> phase_grad.cpt
#  gmt grdimage xphase.grd $scale -Cphase_grad.cpt -X.2i -Y.5i -P > xphase.ps
#  gmt grdimage yphase.grd $scale -Cphase_grad.cpt -X.2i -Y.5i -P > yphase.ps
#
  mv mask.grd tmp.grd 
  gmt grdmath tmp.grd FLIPUD = mask.grd
#
# delete files
 rm tmp.grd tmp2.grd ximag.grd yimag.grd xreal.grd yreal.grd 
