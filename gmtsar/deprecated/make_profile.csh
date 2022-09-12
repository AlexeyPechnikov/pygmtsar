#!/bin/csh -e
#       $Id$
#
#  Kurt Feigl 20150807
#
# make profile across strike
# 
if ($#argv != 4) then
    echo "make profile across strike"
    echo ""
    echo ""
    echo "Usage: $0:t longitude_in_degrees latitude_in_degrees strike_in_degrees_CW_fromN step_size_in_degrees"
    echo ""
    echo ""
    echo ""
    exit 1
endif

# longitude of crossing point
set maxx = $1

# latitude of crossing point
set maxy = $2

# strike in degrees clockwise from north
set dir_strike = $3
set dir_cross  = `gmt math -Q $dir_strike 90. ADD 360. MOD =`

# step size in latitude
#set step = `gmt math -Q 1 360 DIV =`
set step = $4

# draw a line through the maximum (units will be km)
#echo $maxx $maxy
#gmt project -G0.1 -C$maxx/$maxy -A90. -Q -L-50/50  >! profile.rsp

##gmt project -G0.1 -C$maxx/$maxy -A$dir_cross -Q -L-50/50  >! profile.rsp

# make 5 more lines
#set newy = `gmt math -Q 1 360 DIV 2. MUL $maxy EXCH SUB =`
set newxy = "$maxx/$maxy"
#echo newxy is $newxy

\rm -f profile.rsp
touch profile.rsp
set j = 1
while ( $j <= 50 )
    #set newy = `gmt math -Q 1 360 DIV $newy ADD =`
    # calculate coordinates of a point located 0.1 km along strike
    set newxy = `gmt project -G0.1 -C$newxy -A$dir_strike -Q -L0/0.1 | awk 'NR==2{print $1"/"$2}'`
    #echo newxy is $newxy
    echo '>' >> profile.rsp
    gmt project -G0.1 -C$newxy -A$dir_cross -Q -L-50/50 >> profile.rsp
    # calculate coordinates of a point located 1 km along strike
    set newxy = `gmt project -G0.1 -C$newxy -A$dir_strike -Q -L-0.1/0. | awk 'NR==2{print $1"/"$2}'`
    #echo newxy is $newxy
    echo '>' >> profile.rsp
    gmt project -G0.1 -C$newxy -A$dir_cross -Q -L-50/50 >> profile.rsp
    @ j++
end

echo "created:"
\ls -l profile.rsp
exit 0


