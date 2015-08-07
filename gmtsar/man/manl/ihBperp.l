.TH IHBPERP l "06/30/98 \(co 1998 David Sandwell" "ihBperp.l 1.0" "SIOSAR Command"
.SH NAME
ihBperp \- computes a file of perpendicular baselines
.SH SYNOPSIS
.B ihBperp
.I phase.head  baseline.gips
.SH DESCRIPTION
The
.B  ihBperp
command reads the header information from the phase.head file and computes
a matching sized file of perpendicular baselines.  This file can be used
to scale and subtract a reference phase. Here is an example where we
have two files;  
.I phase.gips 
is wrapped phase from an interferogram and 
.I topophase.gips 
is total unwrapped phase for a synthetic interferogram
with a nominal perpendicular baseline of 100m.
.sp
ihBperp phase.gips scale.gips
.br
iha "p = p1 - p2*p3/100." phase.gips topophase.gips scale.gips > diff_phase.gips
.br
.sp 2
.TP 10
.I phase.head
input gips header file (or data file) containing starting and ending baselines as well as
file size information
.TP 10
.I scale.gips
output file of perpendicular baseline in meters
.SH FILES
.ta \w'/opt/gips/tables/cprofile    'u
/opt/gips/tables/cprofile	image header definition file
.SH SEE ALSO
.B ers_baseline(l) plane_wrap(l)
.sp
.SH AUTHOR
David T. Sandwell, Scripps Institution of Oceanography
.SH DIAGNOSTICS
no bugs yet!!
