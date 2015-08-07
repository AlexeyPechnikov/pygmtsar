.TH EXACT_ERS_TIME l "15 JUN 2000, Copyright 2000 David Sandwell" " " "SIOSAR ERS time correction"
.SH NAME
exact_ers_time \- corrects the start time of the first radar echo using the ICU clock
imbedded in the raw data file
.SH SYNOPSIS
.B exact_ers_time
.I SC_clock_start SC_identity icu_start
.SH DESCRIPTION
The
.B exact_ers_time
command reads the approximate start-time and ICU clock counter of the 
first radar echo and computes a new corrected
start time. The approximate time is used to determine the correct PATN file
that, in turn, is used to translate the onboard computer clock counter 
into UTC. The PATN files are generated at ESA from precise tracking data.
.sp
.br
icu - ldr      - exact start time minus app. start time (s)
.br
SC_clock_start - start time of repeat pass
.br
.br
.SH INPUT
.TP 10
.I SC_clock_start
start time of frame in format YYYYDDD.DDDD 
.TP 10
.I SC_identity
spacecraft number of reference (1)-ERS-1 (2)-ERS2
.TP 10
.I icu_start
onboard clock counter of first record in frame
.sp
.SH INPUT
.br
icu - ldr      - exact start time minus start time (s)
.br
SC_clock_start - exact start time of frame
.br
.SH EXAMPLES
Compute exact start time for frames in the Salton Sea area. These are recent ERS-2
data and the largest difference is 2.7 milliseconds.
.br
.sp
exact_ers_time  1999169.7658487617 2 3555329853.933
.br
icu - ldr =  0.000167
.br
SC_clock_start          = 1999169.76584876352
.br
.sp
exact_ers_time  1999204.7658464119 2 34506524.933
.br
icu - ldr = -0.000049
.br
SC_clock_start          = 1999204.76584641146
.br
.sp
exact_ers_time  1999274.7658527894 2 1582794635.714
.br
icu - ldr =  0.002735
.br
SC_clock_start          = 1999274.76585282106
.br
.sp
exact_ers_time  1999344.7658110070 2 3131081524.933
.br
icu - ldr =  0.000059
.br
SC_clock_start          = 1999344.76581100770
.br
.sp
exact_ers_time  2000014.7658063078 2 3905225265.714
.br
icu - ldr =  0.002198
.br
SC_clock_start          = 2000014.76580633316
.br
.sp
exact_ers_time  2000049.7658661806 2 384403147.933
.br
icu - ldr =  0.001426
.br
SC_clock_start          = 2000049.76586619718
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
