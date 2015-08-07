.TH ERS_LINE_FIXER l "08/28/97 \(co JPL" "ers_line_fixer.l 1.0" 
.SH NAME
ers_line_fixer \- fixes missing lines in raw ERS SAR data 
.SH SYNOPSIS
.B  ers_line_fixer
.I -a -h header -l line -s
.I station ifile ofile
.SH DESCRIPTION
The
.B  ers_line_fixer
fixes missing lines in raw SAR data.

.TP 10
.I -a 
align window sampling times
.TP 10
.I -h 
header length (bytes)
.TP 10
.I -l 
line length in bytes
.TP 10
.I -s 
station - either DPAF/ESRIN, CO, EIC, UK, CCRS, or unknown
.TP 10
.I input_file
.TP 10
.I output_file
.sp
.TP 0
.B example: ers_line_fixer -s CCRS input.raw out.fix

.SH SEE ALSO
read_ccrs(l)
CCT.format* files in /opt/siosar/src/pre_proc/read_ccrs/SAR_tape_data_description
