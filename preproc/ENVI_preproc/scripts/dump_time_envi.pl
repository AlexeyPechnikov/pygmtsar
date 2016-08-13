#!/usr/bin/perl

use strict;
#use warnings;

if ($#ARGV != 0) {
# print $#ARGV;  # very strange, this number should be 0 if no input, but here -1
 print "\nusage: dump_time_envi.pl envisat.baq/N1 \n\n";
 print "input:  envisat.raw: envisat raw data in .N1 or .baq \n";
 print "output: print out start time and end time \n\n";
 exit;
}

my %month = (
	"JAN" => "01",
	"FEB" => "02",
	"MAR" => "03", 
	"APR" => "04",
	"MAY" => "05",
	"JUN" => "06",
	"JUL" => "07",
	"AUG" => "08",
	"SEP" => "09", 
	"OCT" => "10", 
	"NOV" => "11", 
	"DEC" => "12");
my @text = ();
my $file = $ARGV[0]; 
my $file2 = "tmp_time"; # print out so other program can do julian day conversion
my $start = "";
my $stop = "";
my $date = "";
my $time = "";
my $year = "";
my $jday = "";
my $ptime ="";
my $ftime0 ="";
my $ftimey ="";
my $ftime2 ="";
my @tmp = ();
my $n = 0;

@text = `head -15 $file`;

LINE:foreach(@text){
    if($_ =~ /sensing_start="(.*)"/i){
	$start = $1;
	$n++;
    }	
    if($_ =~ /sensing_stop="(.*)"/i){
	$stop = $1;
	$n++;
    }
    last LINE if($n >= 2) ;
}
close File;

@tmp = split /\s+/, $start;
$date = $tmp[0];
$time = $tmp[1];
@tmp = split /-/, $date;
$year = $tmp[2];
$ftime2 = $tmp[2]*10000+$month{$tmp[1]}*100+$tmp[0];
$jday = juldate($tmp[2],$month{$tmp[1]},$tmp[0]); 
@tmp = split /:/, $time;
$ptime = ($tmp[0] + $tmp[1]/60 + $tmp[2]/3600)/24;
$ftime0 = $jday+$ptime;
$ftimey = $year*1000+$jday+$ptime;
print "SENSING_START0 ";
print "$ftime0 \n";
print "SENSING_STARTY ";
print "$ftimey \n";
print "START2 ";
print "$ftime2 \n";

@tmp = split /\s+/, $stop;
$date = $tmp[0];
$time = $tmp[1];
@tmp = split /-/, $date;
$year = $tmp[2];
$ftime2 = $tmp[2]*10000+$month{$tmp[1]}*100+$tmp[0];
$jday = juldate($tmp[2],$month{$tmp[1]},$tmp[0]);
@tmp = split /:/, $time;
$ptime = ($tmp[0] + $tmp[1]/60 + $tmp[2]/3600)/24;
$ftime0 = $jday+$ptime;
$ftimey = $year*1000+$jday+$ptime;
print "SENSING_STOP0 ";
print "$ftime0 \n";
print "SENSING_STOPY ";
print "$ftimey \n";
print "STOP2 ";
print "$ftime2 \n";

`rm -f tmp_time`;

sub juldate {
 my $iyear = $_[0];
 my $imon = $_[1];
 my $iday = $_[2];
 my @jsum = (0,31,59,90,120,151,181,212,243,273,304,334);
 my @jsum2= (0,31,60,91,121,152,182,213,244,274,305,335);
 if ($iyear%400 == 0) {      #leap year
     $jday=$jsum2[$imon-1]+$iday;
 } 
 elsif ($iyear%100 == 0) {  #not a leap year
     $jday=$jsum[$imon-1]+$iday;
 }
 elsif ($iyear%4 == 0) {    #leap year
     $jday=$jsum2[$imon-1]+$iday;
 } else  {                  #not a leap year
     $jday=$jsum[$imon-1]+$iday;
 }
 
 return $jday;
}

