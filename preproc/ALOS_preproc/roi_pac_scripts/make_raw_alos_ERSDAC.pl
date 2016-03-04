#!/usr/bin/perl
### make_raw_alos.pl

use Env qw(INT_SCR INT_BIN MY_BIN);
use lib "$INT_SCR";  #### Location of Generic.pm
use Generic;
use POSIX qw(ceil floor);

sub Usage {
  print STDERR <<END;

Usage: make_raw_alos.pl alos_file_root [outname] [fbd2fbs_conversion] [new_ws_code]
  alos_file_root     : prefix of Level 0 ALOS data file(s) (e.g., IMG-HH-ALP)
  outname            : output file name; YYMMDD is default

Function: Creates I<outname>.raw and I<outname>.raw.rsc from imagery files

make_raw_alos_ERSDAC.pl
put together by Rob Mellors (7/2010, SDSU) but based on make_raw_alos.pl
by Yuri Fialko with modifications by Eric Fielding and Zhenlong Li
software inputs from Dave Sandwell, Rob Mellors, and Matt Wei
* does not handle ALOS_merge or FBS2FBD/FBD2FBS
* requires version of ALOS_pre_process that create roi_pac output (-roi)
* note that ERSDAC data already has a .raw suffix

 Example:
make_raw_alos_ERSDAC.pl PASL10C0609250613491002120000

(requires PASL10C0609250613491002120000.raw and PASL10C0609250613491002120000.ldr)
END
  exit 1;
}

@ARGV >= 1  or Usage();
@args = @ARGV;

$alos_file_prefix   = shift;
$outname            = shift ;

# get name of image file
# ERSDAC uses *.raw to identify image files
# need to avoid confusion with outname.raw
# RJM
$imagery = `ls $alos_file_prefix.raw` or die "No Imagery files\n";
$leaderfile = `ls $alos_file_prefix.ldr` or die "No Leaderfile \n";
chop($imagery);
chop($leaderfile);

#################
Message "Checking I/O";
#################

Log ("make_raw_alos.pl", @args);
 Message "ALOS_pre_process $imagery $leaderfile";
 `ALOS_pre_process $imagery $leaderfile -nodopp -roi`;
#
# ALOS_pre_process (with -roi flag as of 7/2010) creates:
# ALOS_FILE_PREFIX.PRM 		- PRM header file
# ALOS_FILE_PREFIX.raw.raw 	- reformatted raw data
# tmp.DATE.raw.rsc		- partial rsc filr for raw data
# hdr_data_points_DATE.rsc 	- orbit hdr file 
#              DATE is yymmdd format
#
# roi_pac does not use the ALOS_pre_process doppler 
# information so run with -nodopp flag (faster)
#
# read data from tmp rsc file name
$tmp_raw_rsc = `ls tmp.*.raw.rsc`;
chop($tmp_raw_rsc);
($tmp,$date,$tmp,$tmp) = split(/\./,$tmp_raw_rsc);
#
# copy rsc from ALOS_pre_process to tmp_IMAGERY.raw.rsc 
`cp $tmp_raw_rsc tmp_IMAGERY.raw.rsc`;
#
# rename reformatted raw file (with .raw.raw suffix for ERSDAC !)
# to tmp_IMAGERY.raw
`\mv $imagery.raw tmp_IMAGERY.raw`;
#
# orbit information written by ALOS_pre_process
$HDR = "hdr_data_points_".$date.".rsc";
#
# read parameters from rsc file created by ALOS_pre_process
#
$day   = Use_rsc  "tmp_IMAGERY.raw read FIRST_LINE_DAY_OF_MONTH";
$month = Use_rsc  "tmp_IMAGERY.raw read FIRST_LINE_MONTH_OF_YEAR";
$year  = Use_rsc  "tmp_IMAGERY.raw read FIRST_LINE_YEAR";
$first_line_utc = Use_rsc "tmp_IMAGERY.raw read FIRST_LINE_UTC";
$center_utc = Use_rsc "tmp_IMAGERY.raw read CENTER_LINE_UTC";
$last_line_utc = Use_rsc "tmp_IMAGERY.raw read LAST_LINE_UTC";
$equatorial_radius = Use_rsc "tmp_IMAGERY.raw read EQUATORIAL_RADIUS";
$file_length = Use_rsc "tmp_IMAGERY.raw read FILE_LENGTH";

$orbit_type = "HDR";
$sat = "ALOS";

$outname = $date;

Message "Using Orbit Information";

($q1,$q2,$Lat,$Lon,$height_mid, $x0, $y0, $z0, $vx0, $vy0,$vz0) = split /\s+/,
    `$INT_SCR/state_vector.pl $year$month$day $center_utc $sat $orbit_type $outname`;
Status "state_vector.pl";

$pi   = atan2(1,1)*4;

if ($orbit_type eq "HDR"){
	&calc_height_GRS80;
	$height_mid=$H;
        }
#
# used original commetns to name
#
&calc_radius_WGS84;

($q1,$q2,$q3,$q4,$height_top, $x0, $y0, $z0, $vx, $vy,$vz) = split /\s+/,
    `$INT_SCR/state_vector.pl $year$month$day $first_line_utc $sat $orbit_type $outname`;
Status "state_vector.pl";

if ($orbit_type eq "HDR"){
	&calc_height_GRS80;
	$height_top=$H;
        }

$height_dt=($height_mid-$height_top)/($center_utc-$first_line_utc);

if ($vz0 > 0) {$orbit_direction =  "ascending";}
else          {$orbit_direction = "descending";}

$velocity_mid = sqrt($vx0**2 + $vy0**2 + $vz0**2);

$Latd = $Lat*180.0 / $pi;
$Lond = $Lon*180.0 / $pi;
$hdgd = $hdg*180.0 / $pi;

print STDERR "vel $vel $velocity_mid \n";

# write it out to the rsc file
Use_rsc "tmp_IMAGERY.raw write HEIGHT_TOP   $height_top";
Use_rsc "tmp_IMAGERY.raw write HEIGHT       $height_mid";
Use_rsc "tmp_IMAGERY.raw write HEIGHT_DT    $height_dt";
Use_rsc "tmp_IMAGERY.raw write VELOCITY     $velocity_mid";
Use_rsc "tmp_IMAGERY.raw write LATITUDE     $Latd";
Use_rsc "tmp_IMAGERY.raw write LONGITUDE    $Lond";
Use_rsc "tmp_IMAGERY.raw write HEADING      $hdgd";
Use_rsc "tmp_IMAGERY.raw write EQUATORIAL_RADIUS   $equatorial_radius";
Use_rsc "tmp_IMAGERY.raw write ECCENTRICITY_SQUARED $e2";
Use_rsc "tmp_IMAGERY.raw write EARTH_EAST_RADIUS $N";
Use_rsc "tmp_IMAGERY.raw write EARTH_NORTH_RADIUS $M";
Use_rsc "tmp_IMAGERY.raw write EARTH_RADIUS $earth_radius_mid";
Use_rsc "tmp_IMAGERY.raw write ORBIT_DIRECTION $orbit_direction";
#
################################
#Message "Doppler Computation";
################################
#
## use external routine EJF
system "$INT_SCR/scan_doppler.pl tmp_IMAGERY $file_length";
#
## move to final file name
`mv tmp_IMAGERY.raw.rsc             ${outname}.raw.rsc`;
`mv tmp_IMAGERY.raw                 ${outname}.raw`;
#
##########################
#Message "Raw data ready for processing";
##########################
#
exit 0;
#--------------------------------------------------------------------
sub calc_height_GRS80{ 
  $ae    = 6378137;             #GRS80 reference ellipsoid
  $flat  = 1/298.257223563;	# looks like WGS84 ? 
  $r     = sqrt($x0**2+$y0**2+$z0**2);
  $r1    = sqrt($x0**2+$y0**2);
  $Lat   = atan2($z0,$r1);
  $Lon   = atan2($y0,$x0);
  $H     = $r-$ae;

  for ($i=1; $i<7; $i++){
    $N      = $ae/(sqrt(1-$flat*(2-$flat)*sin($Lat)**2));
    $TanLat = $z0/$r1/(1-(2-$flat)*$flat*$N/($N+$H));
    $Lat    = atan2($TanLat,1);
    $H      = $r1/cos($Lat)-$N;
  }
# Message "calc_height_GRS80 $H"; 
}
#--------------------------------------------------------------------
sub calc_radius_WGS84 {
   $ae   = 6378137;                        #WGS84 reference ellipsoid
   $flat = 1./298.257223563;
   $N    = $ae/sqrt(1-$flat*(2-$flat)*sin($Lat)**2);
   $re_mid=$N;

   $ve  = -sin($Lon) * $vx0 + cos($Lon) * $vy0;
   $vn  = -sin($Lat) * cos($Lon) * $vx0 - sin($Lat) * sin($Lon) * $vy0 + cos($Lat) * $vz0;
   $hdg = atan2($ve,$vn);
   $e2  = $flat * (2-$flat);
   $M   = $ae * (1-$e2) / (sqrt(1-$e2 * sin($Lat)**2))**3;
   $earth_radius_mid = $N*$M/($N*(cos($hdg))**2+$M*(sin($hdg))**2);

#   Message "calc_radius_WGS84 $earth_radius_mid";
}
#--------------------------------------------------------------------
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
#--------------------------------------------------------------------

=pod

=head1 USAGE

B<make_raw_alos.pl_ERSDAC> I< alos_file_prefix_root >

=head1 FUNCTION

Creates I<date>.raw and I<date>.raw.rsc from ALOS ERSDAC imagery files

=head1 ROUTINES CALLED

ALOS_pre_process

state_vector.pl

=head1 CALLED BY

none

=head1 FILES USED

PAS*.raw

PAS*.ldr

=head1 FILES CREATED

I<date>.raw

I<date>.raw.rsc

I<date>_parse_line.out

shift.out

shift.out.rsc

=head1 HISTORY

Perl  Script : Yuri Fialko 07/11/2007
Use ALOS_pre_processor and ALOS_merge to pre-process and merge a set of ALOS frames; Eric Fielding, on 9 Oct 2007
Use ALOS_fbd2fbs and ALOS_fbs2fbd to convert between FBD and FBS images; Zhenhong Li, on 19 Oct 2007
Butchered by RJM to handle ERSDAC format - removed ALOS_fbd2fbs ALOS_fbs2fbd and ALOS_merge.... 27 July 2010

=head1 LAST UPDATE

2010/07/27 R Mellors SDSU

=cut
