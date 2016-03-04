#!/usr/bin/perl
use strict;

if ($#ARGV < 1 ) {
# print $#ARGV;  # very strange, this number should be 0 if no input, but here -1
 print "\nusage: dump_orbit_ers.pl orbit_frame orbdir\n\n";
 print "input:  orbit_frame.PRM\n";
 print "        orbdir directory with PRC orbit data \n";
 print "output: orbit_frame.LED store the orbit information, used by ERS_baseline \n\n";
 print "        year day_of_year sec_of_day x y z vx vy vz \n\n";
 exit;
}

#===========================================
# possible bug: when sensing is in the last 
# day of the orbit file. It will find incorrect file.
# For example: e2_22526_2925.LED 
#===========================================

my $prm_file = $ARGV[0].".PRM";
my $led_file = $ARGV[0].".LED";
my $orb_dir = $ARGV[1];

print "---------- orbit search for $ARGV[0] started ----------\n";
####################### grep start and stop time from PRM ###############
my $midtime;

my $temp1 = `grep SC_clock_start $prm_file`;
my $temp2 = `grep SC_clock_stop $prm_file`;
my @temp3 = split(/=/,$temp1);
my @temp4 = split(/=/,$temp2);
my $temp5 = ($temp3[1]+$temp4[1])/2;
my $midtime = jul2ymd_PRM($temp5);

#print "$temp3[1], $temp4[1], $temp5, $midtime \n";


####################### Convert time to days since 1.1.2000 12h ############
my $midtime_julday;
my $midtime_ymd;
my $midtime_y11;
my $first_year;
my $first_d2y; # day of the year for the image

if($midtime =~ /(\d*)-(\d*)-(\d*)/){
  $midtime_ymd = $3.$2.$1;
  $midtime_y11 = $3."0101";
  $first_year = $3;
}
print "Midtime: $midtime; YMD $midtime_ymd; Y11 $midtime_y11; \n";

$midtime_julday = ymd2jul($midtime_ymd);
$first_d2y = $midtime_julday - ymd2jul($midtime_y11) + 1;
print "Day of the year: $first_d2y \n";

my $ref_day = 20000101;
my $ref_julday = ymd2jul($ref_day);

my $julday = ($midtime_julday - $ref_julday)*10 - 5; # unit in PRC is 0.1 day from noon 1.1.2000 in TDT

print "days from 1.1.2000 12 h in UTC: $julday/10 \n";

###################### find the right orbit file ###########################
my @orb_list = `cat $orb_dir/arclist`;
my $list_size = @orb_list;
my $orb_file;
my $orb_file2;
my $i;
my $orb_case = 0;

for ($i = 0; $i < $list_size; $i++) {
   if ($orb_list[$i] =~ /(\S*)\s*(-?\d*)\s*(-?\d*)/){
      #print "$julday, $2, $3 \n";
      if (($julday > $2) && ($julday <= $3) && ( $3-$julday > 5)){
         $orb_case = 1; # one orb file, not at edge
         $orb_file = $orb_dir."/".$1;
         print "$ARGV[0] \n  Case 1: $julday, $2, $3; Orbit  file: $orb_file\n";
         last;
         #exit;
      }
      elsif (($julday > $2) && ($julday <= $3)) {
         $orb_case = 2;
         print "$ARGV[0] \n  Case 2: $julday, $2, $3; \n";
         if ($i == $list_size-1) {
            $orb_case = 1;
            $orb_file = $orb_dir."/".$1;
            last;
            print "End of arclist, $orb_file\n";
            #exit;
         }
         my $otmp1 = $1;
         my $otmp2 = $3;
         if ($orb_list[$i+1] =~ /(\S*)\s*(-?\d*)\s*(-?\d*)/) {
            if ($otmp2 == $2) {
               $orb_file = $orb_dir."/".$otmp1;
               $orb_file2 = $orb_dir."/".$1;
               print "Two orbit files, $orb_file $orb_file2\n";
            }
            else {
               $orb_case = 1;
               $orb_file = $orb_dir."/".$otmp1;
               print "Case 2 but no continuous next file.\n";
            }
         }
         last;
         #exit;
      }
      else {
      }
   }
}

if ($orb_case == 0) {
  print "$ARGV[0] \n  Case 0: No orbit file exist!";
  exit;
}


###################### Convert time to seconds ############################
my $time_msec;
if($midtime =~ /(\d*):(\d*):(\d*)/){
  $time_msec = ($1*3600+$2*60+$3)*1e6; # middle time in microseconds
}
#print "msec of the day in UTC sensing middle: $time_msec\n";

###################### define time window start 28 minutes in TDT ############
open(DAT,$orb_file) || die("Could not open file $orb_file!");
my @orb_data = <DAT>;
close(DAT);

# grep TDT-UTC
my $tdtutc;
if($orb_data[0] =~ /^.{177}(.{5})/){
   $tdtutc = $1*1e3; # microseconds
   #print "TDTUTC: $tdtutc ms \n";
}

my @orb_data2;
if ($orb_case == 2) {
   open(DAT2,$orb_file2) || die("Could not open file $orb_file!");
   @orb_data2 = <DAT2>;
   close(DAT2);

   # grep TDT-UTC
   my $tdtutc2;
   if($orb_data2[0] =~ /^.{177}(.{5})/){
      $tdtutc2 = $1*1e3; # microseconds
   }

   if ($tdtutc != $tdtutc2) {
      print "TDTUTC: $tdtutc $tdtutc2 \n";
      print "TDTUTC differs for the two orbit files. Special case1 \n";
      exit;
   }
}

my $time_msec_tdt = $time_msec + $tdtutc;
if($time_msec_tdt > 86400000000) {
   $time_msec_tdt = $time_msec_tdt - 86400*1e6;
   $julday = $julday + 10;
}

if($time_msec < 14*60*1e6){
   $first_d2y = $first_d2y - 1;
}

my $julday_start;
my $time_start;
$time_start = $time_msec_tdt - 14*60*1e6; # 14 minutes before middle time
if($time_start < 0){
   $julday_start = $julday - 10;
   $time_start = $time_start + 86400*1e6;
}
else{
   $julday_start = $julday;
}

print "Time window start in TDT (days to 1.1.2000 12h;ms of the day):$julday_start/10;$time_start\n";

#print "check:$time_start ms = $time_msec ms + $tdtutc ms - 14 minutes - 86400s(0or1)\n";

################### search orb_file for 28 minutes ####################
my @tmp7;
my $tmp7size;
my @tmp8;
my $count = 0;
my @out_orb_data;
LINE:foreach(@orb_data){
   if($_ =~ /STTERR/){
      @tmp7 = split /STTERR/,$_;
      $tmp7size = @tmp7 -2;
      @tmp8 = @tmp7[1 .. $tmp7size];
  }
}

if ($orb_case == 2) {
  my @tmp82;
  LINE:foreach(@orb_data2){
   if($_ =~ /STTERR/){
      @tmp7 = split /STTERR/,$_;
      $tmp7size = @tmp7 -2;
      @tmp82 = @tmp7[1 .. $tmp7size];
   }         
  }
  my @tmp81 = @tmp8;
  @tmp8 = (@tmp81,@tmp82);
}

my @out_orb_data;
my @tmp9;
my $count = 0;
my $first_s;
my $n_sec;
LINE:foreach(@tmp8){
   if($_ =~ /^.{8}(.{6})(.{11})(.{12})(.{12})(.{12})(.{11})(.{11})(.{11})/){ # ask guoguo
      @tmp9 = ($1,$2,$3,$4,$5,$6,$7,$8);
      $tmp9[0] =~ s/^\s*//;
      $tmp9[0] =~ s/\s*$//;
      #print "$tmp9[0] $tmp9[1] $tmp9[2] $tmp9[3] $tmp9[4]\n";
      if( ($tmp9[0] eq $julday_start) && ($tmp9[1] > $time_start) && ($count < 56)){ # 28 minutes, 30s each
            $count = $count+1;
            if($count eq '1'){
               #print "First time in TDT: $julday_start $tmp9[1] ms (check with start time)\n";
               $first_s = ($tmp9[1] - $tdtutc)/1e6; # convert to UTC
               if($first_s < 0){
                  $first_s = $first_s +86400;
               }
            }
            $n_sec = $first_s + 30*($count-1);
            #$n_sec = ($tmp9[1] - $tdtutc)/1e6;
            $tmp9[2] = sprintf("% .3f",$tmp9[2]/1e3); # change unit to meters
            $tmp9[3] = sprintf("% .3f",$tmp9[3]/1e3);
            $tmp9[4] = sprintf("% .3f",$tmp9[4]/1e3);
            $tmp9[5] = sprintf("% .6f",$tmp9[5]/1e6);
            $tmp9[6] = sprintf("% .6f",$tmp9[6]/1e6);
            $tmp9[7] = sprintf("% .6f",$tmp9[7]/1e6);
            push(@out_orb_data,"$first_year $first_d2y $n_sec $tmp9[2] $tmp9[3] $tmp9[4] $tmp9[5] $tmp9[6] $tmp9[7]\n");
      }
   }
}

print "TDTUTC: $tdtutc ms \n";

###################### write to the output file ###########################
unshift(@out_orb_data,"$count $first_year $first_d2y $first_s 30 \n");
open(DAT2,">$led_file") || die("Could not open file!");
print DAT2 @out_orb_data;
close(DAT2);
#print "56 orbit data was printed to $led_file\n";
print "------- Orbit data finished --------\n";

################## subroutine to convert time format #################
sub ymd2jul {
	 # convert year month day to integer Julian day
	 #print "$_[0]\n";
         # input
	 my $iyear = substr($_[0],0,4);
         my $imonth = substr($_[0],4,2);
	 my $iday = substr($_[0],6,2);
         #print "$iyear $imonth $iday \n";

	 # algrithm from wiki
         my $a = int((14-$imonth)/12);
         my $y = $iyear + 4800 - $a;
         my $m = $imonth + 12*$a - 3;
         my $jday2day = $iday + int((153*$m+2)/5) + 365*$y + int($y/4) - int($y/100) + int($y/400) - 32045;

	 return $jday2day;
}

sub jul2ymd_PRM {
        # convert Julian day to year month day hour min sec
        #$ymd.hms = jul2ymd_PRM($julday);
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

	if ( ($imon < 10) && ($iday < 10) ) {
	  $ymd = sprintf("%s", "0".$iday."-0".$imon."-".$iyear." ".$ihou.":".$imin.":".$isec);
        } 
	elsif ($imon < 10) {
          $ymd = sprintf("%s", $iday."-0".$imon."-".$iyear." ".$ihou.":".$imin.":".$isec);
        }
        elsif ($iday < 10) {
          $ymd = sprintf("%s", "0".$iday."-".$imon."-".$iyear." ".$ihou.":".$imin.":".$isec);
	} else {
           $ymd = sprintf("%s", $iday."-".$imon."-".$iyear." ".$ihou.":".$imin.":".$isec);
        }
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
