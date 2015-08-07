.TH IHCONV l "04/25/98 \(co 1998 David Sandwell" "ihconv.l 1.0" "SIOSAR Command"
.SH NAME
ihconv \- convolves a 2-D filter with an array
.SH SYNOPSIS
.B ihconv
.I idec jdec filter_file < input.gips  >_out.gips
.SH DESCRIPTION
The
.B ihconv
convolves a filter with an array by skipping through the rows
and columns by idec and jdec. The complete input file is used
for the convolution but output can be decimated in one step which
saves a lot of computer time with respect to the GIPS command 
.B ihbox.

.TP 10
.I idec
skips rows at an increment of idec
.TP 10
.I jdec
skips columns at an increment of jdec
.TP 10
.I filter_file
ascii file of filter coefficients in free format.  The first line of the file
gives the number_of_columns, number_of_rows in the filter (i.e. 2 integers).  The coefficients
follow starting at the top and moving across the first row.
.TP 10
.I input.gips
input file to be filtered (standart input only).  This file must have a gips header attached.
The elements can be almost any type but complex numbers will be converted to
amplitudes prior to filtering.
.TP 10
.I output.gips
filtered and decimated output file (standard output only).  
.SH FILES
.ta \w'/opt/gips/tables/cprofile    'u
/opt/gips/tables/cprofile	image header definition file
.SH SEE ALSO
.B ihbox(l) phased_grad(l)
.sp
.SH AUTHOR
David T. Sandwell. 
.SH DIAGNOSTICS
