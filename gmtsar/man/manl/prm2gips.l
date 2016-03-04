.TH PRM2GIPS l "11/17/96 \(co 1996 David Sandwell" "prm2gips.l 1.0" "SIOSAR Conversion Command"
.SH NAME
prm2gips \- translate a SIOSAR parameter file into a GIPS file header for an SLC image. 
.SH SYNOPSIS
.B prm2gips 
.I file.PRM file.head [R4]
.SH DESCRIPTION
The
.B prm2gips
command reads a parameter PRM file used for focussing raw SAR signal data into
a single look complex (SLC) image.  Parameters are used to create a header for
the SLC image.
.TP 10
.I file.PRM
parameter file for SAR processor 
.TP 10
.I file.head
gips header file
.TP 10
.I [R4]
creates header for real*4 SLC rather than the default integer*2 SLC
.SH FILES
.ta \w'/opt/gips/tables/cprofile    'u
/usr/gips/tables/cprofile	image header definition file
.SH SEE ALSO
.B esarp(l)
.sp
General Image Processing System -- GIPS(tm) by Peter G. Ford,
MIT center for Space Research.
.br
http://delcano.mit.edu/info/gips
.SH AUTHOR
David T. Sandwell, Scripps Institution of Oceanography
.SH DIAGNOSTICS
no bugs yet!!
