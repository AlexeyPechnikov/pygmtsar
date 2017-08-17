#!/bin/csh -f
#       $Id$
#
unset noclobber
#
# Rewritten by Xiaohua XU for GMTSAR, Feb 2016
#
# Matt WEI Jan 28 2010
# Modified by Xiaopeng Feb 8 2010 
# Modified by Matt Feb 12 2010 - zero.hgt and add 0 to names as N01W072.hgt
# Note this script requires a local copy of the SRTM or other global dem data.
# We have setup a web site http://gmtsar.ucsd.edu to prepare and deliver the file dem.grd
#=======================================================================
#  script to make topography for interferograms 
#  The USGS elevations are height above WGS84 so this is OK.
if ($#argv < 2) then
error_message:
 echo " "
 echo " Usage: make_dem.csh lonW lonE latS latN DEM_type" 
 echo " "
 echo "        DEM_type: 1--SRTM1 2--SRTM3 3--ASTER"
 echo " "
 echo "        notes: 1  has not been tested if the tiles cross equator or E180" 
 echo "               2  custom dem.grd also available from http://topex.ucsd.edu/gmtsar "
 echo " "
 exit 1
endif 
#
# for plotting
# 
  set scale = -JX8i
#
# set local variables
if ( -f ~/.quiet ) then
    set V = ""
else
	set V = "-V"
endif

#
  if ($5 == 1) then
    set demdir = "/palsar/DEM/SRTM1/"
    set pixel = 0.000277777778
    set INC = 1s
  else if ($5 == 2) then 
    set demdir = "/palsar/DEM/SRTM3/World/"
    set pixel = 0.000833333333333
    set INC = 3s
  else if ($5 == 3) then 
    set demdir = "/palsar/DEM/ASTER/World/"
    set pixel = 0.000277777778
    set INC = 1s
  else 
    goto error_message
  endif
  set system = `uname -p` 
#
#========================Mosaic topo data===============================
#
# Get the center of the frame from SAT_baseline
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Estimate bound of required topography 

set lon1 = $1
set lon2 = $2
set lat1 = $3
set lat2 = $4
set bound = $1/$2/$3/$4

echo $bound $lon1 $lon2 $lat1 $lat2

#-----------------------------------------------------------------------
# Mosaic topography data
# first, convert "lon1 lon2 lat1 lat2" into "lat_s lat_e lon_s lon_e" 
# that can be used in hgt format 

@ lon10 = `echo "$lon1" | awk '{printf("%d",$1) }'`
@ lon20 = `echo "$lon2" | awk '{printf("%d",$1) }'`
@ lat10 = `echo "$lat1" | awk '{printf("%d",$1) }'`
@ lat20 = `echo "$lat2" | awk '{printf("%d",$1) }'`
echo $lon10 $lon20 $lat10 $lat20

# correct the lon and lat if in the west or south
if ($lon10 < 0 ) @ lon10 = $lon10 - 1
if ($lon20 < 0 ) @ lon20 = $lon20 - 1
if ($lat10 < 0 ) @ lat10 = $lat10 - 1
if ($lat20 < 0 ) @ lat20 = $lat20 - 1

set lon_control = 0
set lat_control = 0
set lon_tmp = $lon10
set lat_tmp = $lat10

rm *tmp*

while ( $lat_tmp <= $lat20 ) 
  echo $lat_tmp
  while ( $lon_tmp <= $lon20 ) 
    # get the hgt file name right
      if ($lat_tmp >= 0) then
        set NS = "N"
        @ lat_s = $lat_tmp
      else if ($lat_tmp < 0 ) then
        set NS = "S"
        @ lat_s = - $lat_tmp
      endif
#
      if ($lon_tmp < 0) then
        set WE = "W"
        @ lon_s = - $lon_tmp
      else if ($lon_tmp >= 0) then
        set WE = "E"
        @ lon_s = $lon_tmp
      endif
#
      if ($lat_s < 10 && $lon_s >= 100) then
         set demtmp = $NS"0"$lat_s$WE$lon_s".hgt"
         set demtmpgrd = $NS"0"$lat_s$WE$lon_s".grd"
      else if ($lat_s < 10 && $lon_s >= 10) then
         set demtmp = $NS"0"$lat_s$WE"0"$lon_s".hgt"
         set demtmpgrd = $NS"0"$lat_s$WE"0"$lon_s".grd"
      else if ($lat_s < 10) then
         set demtmp = $NS"0"$lat_s$WE"00"$lon_s".hgt"
         set demtmpgrd = $NS"0"$lat_s$WE"00"$lon_s".grd"
      else if ($lon_s >= 100) then
         set demtmp = $NS$lat_s$WE$lon_s".hgt"
         set demtmpgrd = $NS$lat_s$WE$lon_s".grd"
      else if ($lon_s >= 10) then
         set demtmp = $NS$lat_s$WE"0"$lon_s".hgt"
         set demtmpgrd = $NS$lat_s$WE"0"$lon_s".grd"
      else
         set demtmp = $NS$lat_s$WE"00"$lon_s".hgt"
         set demtmpgrd = $NS$lat_s$WE"00"$lon_s".grd"
      endif
      set range = `echo "$lon_tmp $lat_tmp" | awk '{print $1"/"$1+1"/"$2"/"$2+1}'` 
      set demlat = $NS$lat_tmp.grd
      echo $demtmp $demlat $range
      
      if (! -e $demdir$demtmp) then   # file does not exist
         set demtmp = "zero.hgt"
      endif
      
      ln -s $demdir$demtmp .

# bytes swap
      if ($system =~ "sparc" || $system =~ "powerpc") then 
        gmt xyz2grd $demtmp -G$demtmpgrd -I$pixel -R$range  -N-32768 -ZTLh $V 
      else if ($system =~ "i686" || $system =~ "i386") then
        gmt xyz2grd $demtmp -G$demtmpgrd -I$pixel -R$range  -N-32768 -ZTLhw $V
      endif
# paste grd files together along longitude
      if (-e tmplon.grd) then
         gmt grdpaste $demtmpgrd tmplon.grd -Gtmplon.grd
      else
         cp $demtmpgrd tmplon.grd
      endif
      
      rm -f $demtmp
      
      @ lon_tmp = $lon_tmp + 1
      if ($lon_tmp >= 180) @ lon_tmp = $lon_tmp - 360
   end

# paste grd files together along latitude
   if (-e tmplat.grd) then
      gmt grdpaste tmplon.grd tmplat.grd -Gtmplat.grd
      rm -f tmplon.grd
   else
      mv -f tmplon.grd tmplat.grd
   endif

   set lon_control = 0
   set lon_tmp = $lon10
   @ lat_tmp = $lat_tmp + 1 
   rm -f $demlat
end

#mv tmplat.grd dem_ortho.grd
gmt grdcut -R$bound tmplat.grd -Gdem_ortho.grd

#
#  add the egm96 geoid model if available
#
if (-e $demdir/geoid.egm96.grd) then   # file does not exist
   echo "adding the egm96 geoid to the dem"
   gmt grdsample $demdir/geoid.egm96.grd `gmt grdinfo -I- dem_ortho.grd` `gmt grdinfo -I dem_ortho.grd` -fg -Gegm96.grd -T
   gmt grdedit egm96.grd `gmt grdinfo -I- dem_ortho.grd` $V
   gmt grdmath dem_ortho.grd egm96.grd ADD = dem.grd
   rm -f egm96.grd
else
   mv -f dem_ortho.grd dem.grd
endif

#
# plot the dem.grd 
# 
  gmt grd2cpt dem.grd -Cgray -Z $V > dem.cpt 
  gmt grdimage dem.grd $scale -Baf -BWSne -P $V -X0.2i -Cdem.cpt -K > dem.ps 
  gmt psscale -Rdem.grd -J -DJTC+w5i/0.2i+h -Cdem.cpt -Baf -O >> dem.ps
  gmt psconvert -Tf -P -Z dem.ps
  echo "DEM map: dem.pdf"
#
#   make a kml file of shaded relief
#
  gmt grdgradient dem.grd -N.7 -A325 -Gdem_grad.grd $V
  gmt grd2cpt dem_grad.grd -Cgray -Z $V > dem_grad.cpt
  grd2kml.csh dem_grad dem_grad.cpt
#
# clean up
# 
  rm -f tmp.log all.grd tmplat.grd
  rm -f dem.cpt *.bb *.eps
  rm -f N*.grd S*.grd dem_grad.grd
