#!/bin/csh -f 
#       $Id$
# Usage: m2s.csh pixel_size_in_meters llpfile
#
# Convert pixel dimension in meters to dx/dy in arc seconds at mean latitude
set pix = $1	# Input pixel dimension in meters
set llp = $2	# lon lat phase binary float file
# 1. Get w e s n in array
set range = (`gmt gmtinfo $2 -bi3f -C`)
# 2. Get mean latitude
set mlat = `gmt math -Q ${range[3]} ${range[4]} ADD 2 DIV  = `
# 3. Get nearest integer 1/2 arc second for latitude (at least 1")
set dy = `gmt math -Q $pix 111195.079734 DIV 3600 MUL 2 MUL RINT 1 MAX 2 DIV  = `
# 4. Get nearest integer 1/2 arc second for longitude at mean latitude (at least 1")
set dx = `gmt math -Q $pix 111195.079734 DIV $mlat COSD DIV 3600 MUL 2 MUL RINT 1 MAX 2 DIV  = `
# Report two -Idx/dy settings: first for actual grid and 2nd for 10 times larger intervals
set inc1 = "${dx}s/${dy}s"
set dx = `gmt math -Q $dx 10 MUL  = `
set dy = `gmt math -Q $dy 10 MUL  = `
set inc2 = "${dx}s/${dy}s"
echo $inc1 $inc2
