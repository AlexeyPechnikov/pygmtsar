.TH PRE_PROC l "04/01/98 \@ 1996 Rob Mellors, Paul Jamason and David Sandwell" "1.0" "SIOSAR Command"
.SH NAME
pre_proc \- reads raw SAR data files and creates a parameter file and a fixed raw data file
.SH SYNOPSIS
.B pre_proc (pre_proc_dpaf or pre_proc_ccrs)
.I sat_orbit_frame  [swst]
.SH DESCRIPTION
Reads 3 files from tape or disk (DPAF-pre_proc_dpaf or CCRS-pre_proc_ccrs format) and creates a parameter file and 
a fixed raw data file that are fed to the SAR processor (esarp). The script
first checks the local disk for two files that it could get from the tape if
necessary (sat_orbit_frame.ldr and sat_orbit_frame.raw). Then the script runs 
.B ers_line_fixer
to create sat_orbit_frame.fix.  The optional parameter 
.I swst
is the sample window start time for proper near-range alignment of raw
data.  The default value is 878 which is good for the radius of the earth
at the latitude of Southern California.  For other latitudes, this number
will need to be changed.
During this processing, the parameter file (sat_orbit_frame.PRM) is constructed.
 
.TP 10
.I sat_orbit_frame
orbit number and frame number of data tape
.TP 10
.I sat_orbit_frame.PRM
parameter file used for the SAR processor (esarp)
.TP 10
.I sat_orbit_frame.fix
fixed raw signal data to be read by the SAR processor
.SH FILES
/opt/siosar/lib/virgin.PRM
.SH SEE ALSO
.B read_ccrs ers_line_fixer read_sarleader read_data_file ers_baseline
.sp
CCT.format* files in SAR_tape_data_description
.sp
.SH AUTHORS
Rob Mellors, Paul Jamason, and David Sandwell, Scripps Institution of Oceanography
.SH DIAGNOSTICS
