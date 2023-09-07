#!/usr/bin/perl -w
use strict;
use Cwd;
use constant PI    => 4 * atan2(1, 1);
use Math::Complex;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#
#   Copyright (c) 2008-2011  Anna Kelbert - All rights reserved
#
# Reads a WS format data file, creates a data file for ModEM inversion.
#
# Last modified: 2012-04-27
#####################################################################

my ($infile, $outfile);

if (@ARGV < 1){
  print "Usage: ./readWSdata.pl InputFile OutputFile\n";
  exit(1);
}
$infile = $ARGV[0];

if (@ARGV < 2){
  $outfile = 'ModEM_input.dat';
} else {
  $outfile = $ARGV[1];
}

my (%easting, %northing);
my (@easting, @northing) = ();
my (%data,%error,%ermap);
my ($site,$period,$direction);

my (@dat,@err) = ();

open (DATA, "$infile") || die("Could not open file $infile!");
my $line = <DATA>;
my ($nsite,$nper,$ncomp) = split(/\s+/,$line);
my @sites;
for (my $i = 0; $i < $nsite; ++$i){
  $sites[$i] = sprintf("%03d",$i);
}
while(<DATA>){
  my $tmp;
  if (/Station_Location/) {
    ($tmp, $direction) = split(/:\s+/,$_);
    if ($direction =~ /N-S/) {
      while (@northing < $nsite){
	$line = <DATA>; chomp $line;
	$line =~ s/\s+//;
	push (@northing, split(/\s+/,$line));
	#my @values = split(/\s+/,$line); print "$#northing: @values\n";
      }
    } elsif ($direction =~ /E-W/) {
      while (@easting < $nsite){
	$line = <DATA>; chomp $line;
	$line =~ s/\s+//;
	push (@easting, split(/\s+/,$line));
	#my @values = split(/\s+/,$line); print "$#easting: @values\n";
      }
    }
  } elsif (/DATA_Period/) {
    chomp;
    ($tmp, $period) = split(/:\s+/,$_);
    for (my $i = 0; $i < $nsite; ++$i){
      $line = <DATA>; chomp $line;
      $data{$period}{$sites[$i]} = $line;
    }
  } elsif (/ERROR_Period/) {
    chomp;
    ($tmp, $period) = split(/:\s+/,$_);
    for (my $i = 0; $i < $nsite; ++$i){
      $line = <DATA>; chomp $line;
      $error{$period}{$sites[$i]} = $line;
    }
  } elsif (/ERMAP_Period/) {
    chomp;
    ($tmp, $period) = split(/:\s+/,$_);
    for (my $i = 0; $i < $nsite; ++$i){
      $line = <DATA>; chomp $line;
      $ermap{$period}{$sites[$i]} = $line;
    }
  }
}
close(DATA);

#my @origin = (795244, 728840);
my @origin = (0, 0);
my $orient = 0.0;

for (my $i = 0; $i < $nsite; ++$i){
  $easting{$sites[$i]} =  $easting[$i] + $origin[0];
  $northing{$sites[$i]} =  $northing[$i] + $origin[1];
}

my $sites = @sites;
my @periods = sort {$a <=> $b} keys %data;
my $periods = @periods;
print "$periods periods processed\n";

# set the error floor of 5%
foreach $period (@periods) {
	foreach $site (sort @sites) {
		$data{$period}{$site} =~ s/\s+//;
		@dat = split(/\s+/,$data{$period}{$site});
#		$error{$period}{$site} =~ s/\s+//;
#      	@err = split(/\s+/,$error{$period}{$site});
    	my $Zxy = $dat[2] + $dat[3]*i;
    	my $Zyx = $dat[4] + $dat[5]*i;
    	my $err = 0.05 * sqrt(abs($Zxy*$Zyx));
    	for (my $i = 0; $i <= $#dat; ++$i){
      		$err[$i] = $err;
    	}
#       WRONG:
#       for (my $i = 0; $i <= $#dat; ++$i) {
#      		$err[$i] = max($dat[$i]*0.05, $err[$i]);
#      	}
      	$error{$period}{$site} = "@err";
	}
}

my @periodsubset;
my $j = 0;
for (my $i = 0; $i <= $#periods; ++$j){
  $periodsubset[$j] = $periods[$i];
  $i += 1;
}
$periods = @periodsubset;
print "using $periods periods\n";

# write a ModEM data file
# example header:
## Comment line
## Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
#> Full_Impedance
#> exp(-i\omega t)
#> Ohm
#> 0.00
#> 0.000 0.000
#> 2 3
my $comp = "ZXX ZXY ZYX ZYY";
my $lat = 0.00;
my $lon = 0.00;
open (DATA, "> $outfile");
print DATA "# Created from $infile using ./readWSdata.pl using 5% errors\n";
print DATA "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n";
print DATA "> Full_Impedance\n";
print DATA "> exp(-i\\omega t)\n";
print DATA "> Ohm\n";
printf DATA "> %6.2f\n", $orient;
printf DATA "> %7.3f %7.3f\n", @origin;
print DATA "> $periods $sites\n";
foreach $site (@sites){
  foreach $period (@periods){
    my @comp = split(/\s+/,$comp);
    my @dat = split(/\s+/,$data{$period}{$site});
    my @err = split(/\s+/,$error{$period}{$site});
    for (my $i = 0; $i <= $#comp; ++$i) {
      printf DATA "%12.6E ",$period; # transmitter
      printf DATA "%s %8.3f %8.3f ",$site,$lat,$lon; # receiver
      printf DATA "%12.3f %12.3f %12.3f ",$northing{$site},$easting{$site},0.0; # receiver x,y,z
      printf DATA "%s %15.6E %15.6E %15.6E\n",$comp[$i],$dat[2*$i],$dat[2*$i+1],$err[2*$i]; # data
    }
  };
}
close (DATA);
