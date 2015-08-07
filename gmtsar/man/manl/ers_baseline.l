.TH ERS_BASELINE l "05 MAR 2000, Copyright 1997 David Sandwell" " " "SIOSAR Baseline Calculator"
.SH NAME
ers_baseline \- calculates baseline parameters from a precise ERS-1/2 orbit
.SH SYNOPSIS
.B ers_baseline
.I ref_t0 ref_tf ref# rep_t0 rep_tf rep# xt yt zt
.SH DESCRIPTION
The
.B ers_baseline
Command reads start-time and stop-time of reference and repeat ERS-1/2 
frames.
Precise orbit information,
stored in ODR files, is used to calculate the position
of the baseline between the reference and repeat orbits at the start
and end of the frame.
The program returns the following information needed for the SAR
processor as well as for interferometry: 
.sp
.br
earth_radius -  radius of ellipse beneath repeat orbit (m)
.br
sc_height    -  height of satellite above ellipse (m)
.br
sc_velocity  -  equivalent velocity for SAR processor (m/s)
.br
baseline_start-  length of baseline (m) at start of frame
.br
alpha_start  -  baseline alpha (deg) at start of frame
.br
baseline_end -  length of baseline (m) at end of frame
.br
alpha_end    -  baseline alpha (deg) at end of frame
.br
xshift       -  number of range cells to allign rep. to ref.
.br
yshift       -  number of echos to allign rep. to ref.
.br 
(Note that if the reference and repeat times are the same then only
the relevant part of the information is returned.)
.br
.SH INPUT
.TP 10
.I ref_t0
starting time of reference frame in format YYYYDDD.DD 
.TP 10
.I ref_tf
ending time of reference frame in format YYYYDDD.DD
.TP 10
.I ref#
spacecraft number of reference (1)-ERS-1 (2)-ERS2
.TP 10
.I rep_t0
starting time of repeat frame in format YYYYDDD.DD 
.TP 10
.I rep_tf
ending time of repeat frame in format YYYYDDD.DD
.TP 10
.I rep#
spacecraft number of repeat (1)-ERS-1 (2)-ERS2
.sp
.SH EXAMPLE
Compute orbit parameters for an interferogram spanning the Landers 1992.
.br
.sp
%ers_baseline 1992220.7697722 1992220.7699099 1 1992115.7695606 1992115.7696984 1
.br
.sp
ref_clock_start         = 1992220.7697721999
.br
ref_clock_stop          = 1992220.7699098999
.br
ref_identity            =    1
.br
SC_clock_start          = 1992115.7695605999
.br
SC_clock_stop           = 1992115.7696984001
.br
SC_identity             =    1
.br
earth_radius            = 6371413.7408
.br
SC_height               =    788146.89
.br
SC_vel                  =  7124.7584
.br
yshift                  =    61
.br
baseline_start          =   150.588
.br
baseline_end            =   146.081
.br
alpha_start             =  153.530
.br
alpha_end               =  152.028
.br
B_parallel              =  -109.18
.br
B_perpendicular         =  -103.72
.br
xshift                  =    13
.br
.sp
Note the two SAR frames reach the reference point at slightly different times corresponding
to an along-track shift of 61.5 echos.  Also note the baseline 
(B_start = 150.588 m, B_end =146.081,  alpha_start = 153.530deg, alpha_end = 152.028deg ) 
closely match the values published in 
Zebker et al., 1994 (B = 146.1 m, alpha = -152 deg).  
.br
.sp
.SH SEE ALSO
.B getorb
software and precise orbit files were provided by Remko Scharroo 
(remko.scharroo@lr.tudelft.nl).
For more up-to-date information on this package and the DUT orbits, netsurf
to:
.br 
http://www.deos.tudelft.nl/ers/precorbs/
.sp
.br
.SH AUTHORS
David T. Sandwell and Karen Watson, Scripps Institution of Oceanography
.SH DIAGNOSTICS
11 NOV 1997 - bug fixed on computation of baselines
.br
03 MAR 2000 - fixed Y2K bugs and expanded year to YYYY
