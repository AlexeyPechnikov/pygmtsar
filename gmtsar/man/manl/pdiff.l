.TH PDIFF l "06/26/97 \(co 1996 David Sandwell" "pdiff.l 1.0" "SIOSAR Command"
.SH NAME
pdiff \- calculates interferogram and subtracts phase due to a curved earth, topography,
and a phase model.
.SH SYNOPSIS
.B pdiff
.I ref.gips rep.gips real.gips imag.gips [topophase.gips] [modelphase.gips]
.SH DESCRIPTION
The
.B pdiff
command reads two matched SLC SAR images, (reference and repeat) and forms
an interferogram.  Two baseline estimates at the start and end of the frame, 
stored in the header of the repeat
image, are used to correct the phase for the "earth flattening".
As an option, a file topophase.gips can be supplied to remove the phase
due to the topography above the local spherical radius of the Earth.
This file must have dimensions that are integral
factors of the SLC images.  For example if the SLC image has 5120 columns
and 28000 rows, topophase.gips could have dimensions 2650 by 7000. The
topographic phase must correspond to a perpendicular baseline of 100 m.
In addition to a topophase, a model-phase can be removed from the interferogram.
This phase is not scaled by the baseline so it represents either ground
deformation or a phase correction due to atmospheric or ionospheric delay.
The dimensions of the model-phase file must match the topophase file.

The real and imaginary components of the 
interferogram are written to separate files. The reference and repeat input files
do not have to have the gips headers attached to the SLC image.  One can
supply just the header files to pdiff as long as the header contains a pointer
to the SLC image.  An example pointer is: 
.br
input_file = 19369_585.SLC 
.TP 10
.I ref.gips
reference SLC image or header of reference SLC image containing pointer to image
.TP 10
.I repeat.gips
repeat SLC image or header of repeat SLC image containing pointer to image
.TP 10
.I real.gips
real component of interferogram
.TP 10
.I imag.gips
imaginary component of interferogram
.TP 10
.I [topophase.gips]
optional file of phase for topography with a 100 m standard baseline
.TP 10
.I [modelphase.gips]
optional file of phase for direct removal from interferogram
.SH FILES
.ta \w'/opt/gips/tables/cprofile    'u
/opt/gips/tables/cprofile	image header definition file
.SH SEE ALSO
.B esarp(l) prm2gips(l) ers_baseline(l) ihBperp(l)
.sp
.SH AUTHOR
Evelyn Price and David T. Sandwell, Scripps Institution of Oceanography
.SH DIAGNOSTICS
no bugs yet!!
