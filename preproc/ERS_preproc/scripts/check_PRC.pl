#!/usr/bin/perl
use strict;

if ($#ARGV < 0 ) {
 print "\nusage: check_PRC.pl PRC_file\n\n";
 print "output: First and last few records\n";
 exit;
}

my $prc_file = $ARGV[0];
my @prc_data;

open(DAT,$prc_file) || die("Could not open file $prc_file!");
my @prc_data = <DAT>;
close(DAT);

my @tmp7;
my $tmp7size;
my @tmp8;
LINE:foreach(@prc_data){
   if($_ =~ /STTERR/){
      @tmp7 = split /STTERR/,$_;
      $tmp7size = @tmp7 -2;
      @tmp8 = @tmp7[1 .. $tmp7size];
  }
}

my @tmp9;
my $count = 0;
my $rsize = @tmp8;

LINE:foreach(@tmp8){
   if($_ =~ /^.{8}(.{6})(.{11})(.{12})(.{12})(.{12})(.{11})(.{11})(.{11})/){ # ask guoguo
      $count = $count + 1;
      @tmp9 = ($1,$2,$3,$4,$5,$6,$7,$8);
      $tmp9[0] =~ s/^\s*//;
      $tmp9[0] =~ s/\s*$//;
      if (($count < 5) || ($count > $rsize -5)) {
          print "$count: $tmp9[0] $tmp9[1] $tmp9[2] $tmp9[3] $tmp9[4]\n";
      }
   } 
}
