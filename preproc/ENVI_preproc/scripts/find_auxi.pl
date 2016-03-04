#!/usr/bin/perl
use strict;

my $auxi_file;
my $time_file;
my @auxi_list;
my @data_time;
my $n = 0;

$auxi_file = $ARGV[0];
$time_file = $ARGV[1];

open(DAT,$auxi_file) || die("Could not open file!");
chomp(@auxi_list = <DAT>);
close(DAT);

open(DAT2,$time_file) || die("Could not open file!");
chomp(@data_time = <DAT2>);
close(DAT2);

LINE:foreach(@auxi_list){
  if($_ =~ /IEC(\d*)_(\d*)_(\d*)_(\d*)_(\d*)_(\d*)/){
      if(($3 < $data_time[0]) && ( $5 > $data_time[1]) ) {
         print $_;
         $n++;
      }
  }
  last LINE if($n >=1);
}
