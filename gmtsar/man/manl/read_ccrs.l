.TH READ_CCRS l "08/28/97 \@ 1996 David Sandwell and Rob Mellors" "1.0" "SIOSAR Command"
.SH NAME
read_ccrs  \- reads files from CCRS SAR data tape and writes to disk 
.SH SYNOPSIS
.B read_ccrs
.I input_device orbit_frame
.SH DESCRIPTION
Reads files from a raw SAR data tape (CCRS format) and writes 3 files to disk.
 
.TP 10
.I input device
tape drive name (example /dev/rmt/0n)
.TP 10
.I vdf 
SAR volume descriptor file - describes what files are on the tape.
.TP 10
.I leader
Mixed ASCII and binary file which contains parameters about the
acquisition of the SAR data. Time information is not included. 
.TP 10
.I data 
SAR data (also has a mixed ASCII and binary header). Time information
is encoded.
.TP 10
.SH FILES
.ta \w'/opt/gips/tables/cprofile    'u
/opt/gips/tables/cprofile	image header definition file
.SH SEE ALSO
.B read_SAR_tape(l)
.sp
CCT.format* files in SAR_tape_data_description
.sp
.SH AUTHORS
David Sandwell and Rob Mellors, Scripps Institution of Oceanography
.SH DIAGNOSTICS
no bugs, just features
