.TH STACK_TABLE l "4/22/99 \(co 1999 David Sandwell" "stack_table.l 1.0" "SIOSAR SAR Processor"
.SH NAME
stack_table \- prepares a table of Bperp vs days into mission. 
.SH SYNOPSIS
.B stack_table
.I master.PRM slave.PRM [GMT]
.SH DESCRIPTION
The
.B stack_table
command reads master and slave parameter files and creates a table
of 
.I sat_orb slave_time slave_days(1992) Bpl Bperp xshift yshift.
The optional GMT parameter reformats the table for use with
the GMT command pstext.  This command is most useful for assembly
of many repeat sar images into a stack for topographic recovery
and change detection.
.TP 10
.I master.PRM
parameter file for master SAR image for the entire stack
.TP 10
.I slave.PRM
parameter file for one of the slave SAR images
.TP 10
.I [GMT]
option to create a textfile suitable for plotting with pstext
.SH EXAMPLE
.br
stack_table 23390_2925.PRM 15240_2925.PRM GMT > temp.txt
.br
stack_table 23390_2925.PRM 15741_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 20384_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 20885_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 21386_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 22388_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 22889_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 23390_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 23891_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 24893_2925.PRM GMT >> temp.txt
.br
stack_table 23390_2925.PRM 25394_2925.PRM GMT >> temp.txt
.br
# 
.br
#  now sort the table
.br
#
.br
sort -n +0 -1 < temp.txt > gmttable.txt
.br
#
.br
pstext gmttable.txt -Jx.004/.0045 -R50/2450/-750/770 -B200/100SWen > plot.ps

.SH SEE ALSO
.B pre_proc_ccrs(l) ers_baseline(l)
.SH AUTHOR
David Sandwell, Scripps Institution of Oceanography. UNIX shell script
.SH DIAGNOSTICS
no bugs yet!!
