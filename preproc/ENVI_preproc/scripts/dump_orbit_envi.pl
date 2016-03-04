#!/usr/bin/perl
use strict;

if ($#ARGV < 3 ) {
# print $#ARGV;  # very strange, this number should be 0 if no input, but here -1
 print "\nusage: dump_orbit_envi.pl start_time stop_time envi.LED orbdir\n\n";
 print "input:  start and stop time in Julian day format\n";
 print "        orbdir directory with Doris orbit data \n";
 print "output: envi.LED store the orbit information, used by ENVI_baseline \n\n";
 print "        year day_of_year sec_of_day x y z vx vy vz \n\n";
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

my $orb_dir = $ARGV[3];

my @orb_list = `ls $orb_dir`;
my @orb_data;
my @out_orb_data;
my $start_time_ymd;
my $stop_time_ymd;
my $tmp_start;
my $tmp_stop;
my $tmp_time;
my $tmp_time_jul;
my $tmp_jul_day;
my $tmp_sec_of_day;
my $start_time_jul;
my $stop_time_jul;
my $start_orb_time;
my $stop_orb_time;
my $orb_file;
my $out_file;
my $n;
my $count;
my $first_time;
my $test_start_time_jul;
my $test_stop_time_jul;
####################### grep start and stop time from PRM ###############
$start_time_jul = $ARGV[0];
$stop_time_jul = $ARGV[1];
$out_file = $ARGV[2];
$start_time_ymd = jul2ymd($start_time_jul);
$stop_time_ymd = jul2ymd($stop_time_jul);
$test_start_time_jul = ymd2jul($start_time_ymd);
$test_stop_time_jul = ymd2jul($stop_time_ymd);
$start_orb_time = ($start_time_jul+$stop_time_jul)/2-14/60/24; # 28 points
$stop_orb_time = ($start_time_jul+$stop_time_jul)/2+14/60/24;

#print "Julian input: $start_time_jul $stop_time_jul \nYMD converted:$start_time_ymd $stop_time_ymd \nJulian converted back:$test_start_time_jul $test_stop_time_jul\n";
#print "orb time range: $start_orb_time $stop_orb_time \n";

####################### Orbit file ######################################

$n = 0;
LINE:foreach(@orb_list){
  if($_ =~ /AXVF-P(\d*)_(\d*)_(\d*)_(\d*)_(\d*)_(\d*)/){
	$tmp_start = $3.".$4";
	$tmp_stop = $5.".$6";
	#print "$tmp_start,$tmp_stop \n";
        if(($tmp_start < $start_time_ymd) && ( $tmp_stop > $stop_time_ymd) ) {
          print "Orbit data: \n".$orb_dir."/".$_;
	  $orb_file = $orb_dir."/".$_;
          $n++;
       }
  }
  last LINE if($n >=1);
}

####################### Read orbit file and print out 28 points #######
open(DAT,$orb_file) || die("Could not open file!");
@orb_data = <DAT>;
close(DAT);

$count = 0;
LINE:foreach(@orb_data){
   if($_ =~ /^(\d*)-(\w*)-(\d*) (\d*):(\d*):(\d*).(\d*) (\S*) (\S*) (\S*) (\S*) (\S*) (\S*) (\S*) (\S*)/){
	$tmp_time = $3.$month{$2}.$1.".".$4.$5.$6;
	$tmp_time_jul = ymd2jul($tmp_time);
	$tmp_jul_day = int($tmp_time_jul)-int($tmp_time_jul/1000)*1000;
	$tmp_sec_of_day = $4*60*60+$5*60+$6;
	#print "$tmp_time_jul\n";
	if (($tmp_time_jul > $start_orb_time) && ($tmp_time_jul < $stop_orb_time)){
	   $count++;
	   if ($count == 1) {
	       $first_time = "$3 $tmp_jul_day $tmp_sec_of_day";
	       push(@out_orb_data,"$3 $tmp_jul_day $tmp_sec_of_day $10 $11 $12 $13 $14 $15\n");
	   }
	   else {
		push(@out_orb_data,"$3 $tmp_jul_day $tmp_sec_of_day $10 $11 $12 $13 $14 $15\n");
	   }
       	}
   }
}

unshift(@out_orb_data,"$count $first_time 60 \n");
open(DAT2,">$out_file") || die("Could not open file!");
print DAT2 @out_orb_data;
close(DAT2);

################## subroutine to convert time format #################
sub ymd2jul {
	 # convert year month day hour min sec to Julian day
	 #$jul.day = ymd2jul($ymd.hms);
	 #print "$_[0]\n";
	 my $iyear = substr($_[0],0,4);
	 my $imon = substr($_[0],4,2);
	 my $iday = substr($_[0],6,2);
	 my $ihou = substr($_[0],9,2);
	 my $imin = substr($_[0],11,2);
	 my $isec = substr($_[0],13,2);
	 my $jday;
	 my $julday;
	 #print "$_ \n$iyear $imon $iday $ihou $imin $isec \n";
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

	 $julday = $iyear*1000 + $jday + $ihou/24 + $imin/60/24 + $isec/60/60/24;
	 return $julday;
}

sub jul2ymd {
	# convert Julian day to year month day hour min sec
	#$ymd.hms = jul2ymd($julday);
	my $iyear = substr($_[0],0,4);
	my $jday = substr($_[0],4,3);
	my $decimal = sprintf("%.8f", ($_[0] - int($_[0]))); 
        my $imon;
        my $iday;
	my $ihou;
	my $imin;
	my $isec;
	my $ymd;
        my @jsum = (0,31,59,90,120,151,181,212,243,273,304,334);
        my @jsum2= (0,31,60,91,121,152,182,213,244,274,305,335);
        if ($iyear%400 == 0) {      #leap year
		$imon = findmonth($jday,@jsum2);
		$iday = $jday-$jsum2[$imon-1];
        }       
        elsif ($iyear%100 == 0) {  #not a leap year
		$imon = findmonth($jday,@jsum);
		$iday = $jday-$jsum[$imon-1];
        }
        elsif ($iyear%4 == 0) {    #leap year
		$imon = findmonth($jday,@jsum2);
		$iday = $jday-$jsum2[$imon-1];
        } else  {                  #not a leap year
		$imon = findmonth($jday,@jsum);
		$iday = $jday-$jsum[$imon-1];
        }

	$ihou = int(24*$decimal);
	$imin = int(60*24*$decimal-60*$ihou);
	$isec = int(60*60*24*$decimal-60*60*$ihou-60*$imin);
	#print "$_ \n$iyear $imon $iday $ihou $imin $isec \n";
	
	$ymd = sprintf("%.6f", $iyear*10000+$imon*100+$iday+$ihou/100+$imin/10000+$isec/1e6);
        return $ymd;
}

sub findmonth {
	#$imon = findmonth($jday,@jsum);
	my $jday = $_[0];
	my @jsum = ($_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7],$_[8],$_[9],$_[10],$_[11],$_[12]);
	#print "$jday @jsum \n";
	my $n = 0;
	my $count = 0;
	my $imonth;
	my $diff;
	LINE:foreach(@jsum){
		$count++;
		$diff = $_ - $jday;
		#print "$count $diff \n";
		if ($diff >= 0){
		   $n++;
		}
	last LINE if($n >=1);
	}
	if ( ($count == 12) & ($diff <= 0) ){
		$imonth = 12;
	} else {
		$imonth = $count-1;
	}
	#print "$n $count $imonth \n";
	return $imonth;
}	
