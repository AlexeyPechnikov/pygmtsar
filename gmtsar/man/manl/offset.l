.TH OFFSET l "04/10/98 \(co 1998 David Sandwell" "offset.l 1.0" "SIOSAR Command"
.SH NAME
offset \- calculates offsets in range and azimuth between master and slave SCL images
.SH SYNOPSIS
.B offset
.I master.head  slave.head offset.dat [-V]
.SH DESCRIPTION
The
.B offset
command reads two SLC SAR images, (master and slave) and computes
row and column offsets needed to align the slave to the master. These offsets
are used by a program 
.B fitoffset 
to generate parameters for the SAR processor 
.B
esarp. Here is an example run:
.sp
.br
offset 22876_585.head 19369_585.head offset.dat
.br
.sp
The second header file (created with the 
.B prm2gips 
command) contains
initial offsets needed to bring the alignment to within 40 pixels.  These
initial offsets can be estimated using 
.B ers_baseline.
.sp
.br 
  xshift = 0
  yshift = 653
.sp
.br
Then one runs 
.B fitoffset 
to refine these parameters.
.sp
.br
fitoffset offset.dat 50
.br
  stretch_r              =    -9.25792E-06
.br
  stretch_a              =     4.23417E-06
.br
  yshift                   =   653
.br
  sub_int_a              =    0.333923
.br
  xshift                   =   -3
.br
  sub_int_r              =    0.806506
.br
.sp

.TP 10
.I master.head
gips header for master image must contain pointer to SLC file. Example pointer is
.br
input_file = 22876_585.SLC
.TP 10
.I slave.head
gips header for slave image must also contain pointer to SLC file.
.TP 10
.I offset.dat
ascii output file containing range and azimuthal offsets for numerous points
in the master image.  SNR of crosscorrelation is also output in the last column.
.TP 10
.I -V
optional verbose mode 
.SH FILES
.ta \w'/opt/gips/tables/cprofile    'u
/opt/gips/tables/cprofile	image header definition file
.SH SEE ALSO
.B esarp(l) prm2gips(l) ers_baseline(1) fitoffset(l)
.sp
.SH AUTHOR
Howard Zebker. Modified by David T. Sandwell. 
.SH DIAGNOSTICS
no bugs yet!!
