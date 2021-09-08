#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# See page 24 "Running Short BAseline Subset (SBAS) Analysis" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# use after GMTSAR.unwrap_orbit.sh and ... next
# ./GMTSAR.sbas_orbit.sh /mnt/GMTSAR asc
set -e

workdir="$1"
orbit="$2"
opts="$3"

cd "$workdir"
cd "$orbit"

rm -fr SBAS/*
cd SBAS

# link the required merged files
ln -f -s ../merge/intf.in            .
ln -f -s ../merge/baseline_table.dat .
ln -f -s ../merge/batch_tops.config  .
ln -f -s ../merge/trans.dat          .
ln -f -s ../merge/*.PRM              .
ln -f -s ../merge/*.LED              .

# Prepare the Input Files
prep_sbas.csh intf.in baseline_table.dat ../merge unwrap.grd corr.grd

# Run SBAS
N=$(wc -l intf.in   | cut -d ' ' -f1)
S=$(wc -l scene.tab | cut -d ' ' -f1)

lon0=$(cat ../reframed/pins.ll | cut -d ' ' -f1 | paste -s -d+ | sed -E 's|(.*)|scale=6; (\1)/2|' | bc -l)
lat0=$(cat ../reframed/pins.ll | cut -d ' ' -f2 | paste -s -d+ | sed -E 's|(.*)|scale=6; (\1)/2|' | bc -l)
elevation=$(gdallocationinfo -valonly -geoloc ../topo/dem.grd "$lon0" "$lat0")

# SAT_look works in local directory files only
# and we have prepared here only one master PRM
satlook=$(echo "$lon0 $lat0 $elevation" | SAT_look *.PRM)
look_E=$(echo "$satlook" | cut -d ' ' -f4)
look_N=$(echo "$satlook" | cut -d ' ' -f5)
look_U=$(echo "$satlook" | cut -d ' ' -f6)
# sbas argument calculation
# 4*a(1) is PI value
incidence=$(echo "scale=3;a(sqrt(${look_E}^2 + ${look_N}^2)/${look_U})*180/(4*a(1))" | bc -l)

# n_columns is the value for x
xdim=$(find -L ../merge -name 'unwrap.grd' | head -n 1 | xargs -I {} -n 1 gmt grdinfo {} | grep n_columns | sed -E 's/(.*) ([0123456789]+)$/\2/')
# n_rows is the value for y
ydim=$(find -L ../merge -name 'unwrap.grd' | head -n 1 | xargs -I {} -n 1 gmt grdinfo {} | grep n_rows    | sed -E 's/(.*) ([0123456789]+)$/\2/')
# x_min
xmin=$(find -L ../merge -name 'unwrap.grd' | head -n 1 | xargs -I {} -n 1 gmt grdinfo {} | grep x_min | sed -E 's/(.*) x_min: ([0123456789]+) (.*)/\2/')
# x_max
xmax=$(find -L ../merge -name 'unwrap.grd' | head -n 1 | xargs -I {} -n 1 gmt grdinfo {} | grep x_max | sed -E 's/(.*) x_max: ([0123456789]+) (.*)/\2/')
# near_range
near_range=$(find -L . -name *.PRM | head -n 1 | xargs -I {} -n 1 grep near_range {} | sed -E 's/(.*) ([0123456789]+\.[0123456789]+)[ \t]*$/\2/')
# sbas argument calculation
rng_samp_rate=$(find -L . -name *.PRM | head -n 1 | xargs -I {} -n 1 grep rng_samp_rate {} | sed -E 's/(.*) ([0123456789]+\.[0123456789]+)[ \t]*$/\2/' | head -n 1)
# radar wavelength
wavelength=$(find -L . -name *.PRM | head -n 1 | xargs -I {} -n 1 grep radar_wavelength {} | sed -E 's/(.*)(0\.[0123456789]+)[ \t]*$/\2/')

# calculation below requires bc utility
rng_pixel_size=$(echo "300000000 / $rng_samp_rate / 2" | bc -l)
rng=$(echo "$rng_pixel_size * ($xmin+$xmax) /2 + $near_range" | bc -l)
# round value
rng=$(printf %.0f $rng)

echo "
DEBUG output to check SBAS arguments:
N=$N
S=$S

lon0=$lon0
lat0=$lat0
elevation=$elevation

look_E=$look_E
look_N=$look_N
look_U=$look_U
incidence=$incidence

xdim=$xdim
ydim=$ydim
xmin=$xmin
xmax=$xmax
near_range=$near_range
rng_samp_rate=$rng_samp_rate
wavelength=$wavelength

rng_pixel_size=$rng_pixel_size
rng=$rng
"

cmd="sbas intf.tab scene.tab "$N" "$S" "$xdim" "$ydim" -range "$rng" -incidence "$incidence" -wavelength "$wavelength" -rms -dem ${opts}"
echo "$cmd"

# notify user via Telegram
telegram_sendmessage.sh "$cmd"

# run the SBAS processing
sbas intf.tab scene.tab "$N" "$S" "$xdim" "$ydim" -range "$rng" -incidence "$incidence" -wavelength "$wavelength" -rms -dem ${opts}

# Post-SBAS Results: Projecting into Latitude/Longitude
proj_ra2ll.csh trans.dat vel.grd vel_ll.grd
gmt grd2cpt vel_ll.grd -T= -Z -Cjet > vel_ll.cpt
grd2kml.csh vel_ll vel_ll.cpt

# attention: proj_ra2ll.csh does not work in parallel mode due to temporary files with the same names
find . -name 'disp_???????.grd' | sed -E 's|(.*).grd|\1.grd \1_ll.grd|' | xargs -n 2 -P 1 proj_ra2ll.csh trans.dat

# define values range of velocity map
range=$(gmt grdinfo vel_ll.grd | grep z_min | sed -E 's/.*z_min: ([-\.0-9]{6}).*z_max: ([-\.0-9]{6}).*/\[\1,\2\]/')
# notify user via Telegram
telegram_sendmediagroup.sh "Velocity Map, ${range} mm/year" vel_ll.png vel_ll.kml
