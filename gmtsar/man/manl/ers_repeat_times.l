.TH ERS_REP_TIMES l "05 MAR 2000, Copyright 1997 David Sandwell" " " "SIOSAR Repeat Frame Alignment"
.SH NAME
ers_repeat_times \- calculates repeat_times of an ERS-2 pass given start and end times
of an ERS-1/2 frame
.SH SYNOPSIS
.B ers_repeat_times
.I ref_t0 ref_tf ref# rep_t0 rep# 
.SH DESCRIPTION
The
.B ers_repeat_times
Command reads start-time and stop-time of reference frame and computes the
corresponding times for a repeat pass. 
Precise orbit information,
stored in ODR files, is used to calculate the position
of the repeat_times corresponding to the times at the start and end of the 
reference orbit.
The program returns the following information needed for makedpaf which
extracts subsets of data from long data takes.
.sp
.br
SC_clock_start - start time of repeat pass
.br
SC_clock_stop  - stop time of repeat pass
.br
SC_ICU_start   - start value of ICU onboard clock counter
.br
SC_ICU_stop    - stop value of ICU onboard clock counter
.br
B_parallel     - approximate parallel baseline
.br
B_perpendicular- approximate perpendicular baseline
.br
xshift       -  number of range cells to allign rep. to ref.
.br 
.br
.SH INPUT
.TP 10
.I ref_t0
starting time of reference frame in format YYYYDDD.DDDD 
.TP 10
.I ref_tf
ending time of reference frame in format YYYYDDD.DDDD
.TP 10
.I ref#
spacecraft number of reference (1)-ERS-1 (2)-ERS2
.TP 10
.I rep_t0
starting day of repeat frame in format YYYYDDD.00 
.TP 10
.I rep#
spacecraft number of repeat (2)-ERS2
.sp
.SH EXAMPLE
Compute repeat times and ICU clock counter for frames in the Salton Sea area.
.br
.sp
ers_repeat_times 1996004.7659481482 1996004.7661548401 1 99169 2
.br
SC_clock_start            = 1999169.7658588917
.br
SC_clock_stop             = 1999169.7660655845
.br
SC_MJD_start            = 18065.7658588917
.br
SC_MJD_stop             = 18065.7660655846
.br
SC_ICU_start             =  3555330077.95
.br
SC_ICU_stop              =  3555334649.67
.br
B_parallel                =  -162.02
.br
B_perpendicular           =  -540.64
.br
xshift                    =    20
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
.SH AUTHOR
David T. Sandwell, Scripps Institution of Oceanography
.SH DIAGNOSTICS
11 NOV 1997 - bug fixed on computation of repeat_times
05 MAR 2000 - fixed Y2K bugs and expanded year to YYYY
