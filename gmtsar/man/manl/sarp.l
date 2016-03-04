.TH ESARP l "4/10/98 \(co 1996 Evelyn Price" "esarp.l 1.0" "SIOSAR SAR Processor"
.SH NAME
esarp \- focusses raw signal data into a single-look complex image (SLC). 
.SH SYNOPSIS
.B esarp 
.I file.PRM file.SLC [R4]
.SH DESCRIPTION
The
.B esarp
command reads a parameter PRM file used provide all of the information needed
tp process the raw signal data data into
a single look complex (SLC) image. Program looks for file.fix, a fixed raw
signal file, for an input stream.
.TP 10
.I file.PRM
parameter file for SAR processor 
.TP 10
.I file.SLC
single-look complex image, 2-integer*2 per pixel
.TP 10
.I [R4]
option to print output as 2-real*4 per pixel
.SH SEE ALSO
.B pre_proc(l) ers_baseline(l) offset(l)
.SH AUTHOR
Evelyn Price and David Sandwell, Scripps Institution of Oceanography. C-code
based on JPL-heritage SAR processor soi.
.SH DIAGNOSTICS
no bugs yet!!
