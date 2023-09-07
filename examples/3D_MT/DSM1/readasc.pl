#!/usr/bin/perl -w
use strict;
use Cwd;
use constant PI    => 4 * atan2(1, 1);
use Math::Complex;

#
#   Copyright (c) 2008  Anna Kelbert - All rights reserved
#
# Reads a bunch of asc input files, creates a single data file for ModEM inversion.
#
# Last modified: 2008-10-16
#####################################################################

my $dir = 'original';
my $file;
my @sites;
opendir(LS, "$dir/") || die "Unable to open the requested directory $dir";
while ($file = readdir(LS)) {
  if ($file =~ /^(\w+).asc$/) {
    push (@sites, $1);
  }
}
closedir(LS);

# read *.asc files
my ($comment1,$comment2);
my (%easting, %northing);
my (%data,%error);
my $site;
foreach $site (sort(@sites)) {
  open (DATA, "$dir/$site.asc");
  my @lines = <DATA>;
  ($comment1,$comment2,@lines) = @lines;
  my @fields = split(/\s+/,$comment1);
  $easting{$site} = $fields[4]; # + 30000;
  $northing{$site} = $fields[7]; # + 30000;
  print "$site $easting{$site} $northing{$site}\n";
  my ($line,@temp);
  foreach $line (@lines){
    my ($period,@data) = split(/\s+/,$line);
    # conversion from Ohm units to (V/m)/T for the Fortran program
    #for (my $i = 0; $i <= $#data; ++$i){
    #  $data[$i] = 10000/(4*PI) * 1000 * $data[$i];
    #}
    @temp = map {sprintf("%.6E",$_)} @data; # convert real array to string array
    $data{$period}{$site} = "@temp"; # join string array values together
    my @error;
    my $value;
    # 0.05*sqrt(abs(Zxy*Zyx))
    my $Zxy = $data[2] + $data[3]*i;
    my $Zyx = $data[4] + $data[5]*i;
    my $err = 0.05 * sqrt(abs($Zxy*$Zyx));
    for (my $i = 0; $i <= $#data; ++$i){
      $error[$i] = $err;
    }
    @temp = map {sprintf("%.6E",$_)} @error;
    $error{$period}{$site} = "@temp";
  }
  close (DATA);
};
my $sites = @sites;
my @periods = sort {$a <=> $b} keys %data;
my $periods = @periods;
print "$periods periods processed\n";

my @periodsubset;
my $j = 0;
for (my $i = 0; $i <= $#periods; ++$j){
  $periodsubset[$j] = $periods[$i];
  $i += 2;
}
$periods = @periodsubset;
print "using $periods periods\n";

my @origin = (0, 0);
my $orient = 0.0;

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
my $period;
open (DATA, "> DSM1.dat");
print DATA "# Synthetic data set DSM1 from Dublin Institute (2008)\n";
print DATA "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n";
print DATA "> Full_Impedance\n";
print DATA "> exp(+i\\omega t)\n";
print DATA "> Ohm\n";
printf DATA "> %6.2f\n", $orient;
printf DATA "> %7.3f %7.3f\n", @origin;
print DATA "> $periods $sites\n";
foreach $site (@sites){
  foreach $period (@periodsubset){
    my @comp = split(/\s+/,$comp);
    my @dat = split(/\s+/,$data{$period}{$site});
    my @err = split(/\s+/,$error{$period}{$site});
    for (my $i = 0; $i <= $#comp; ++$i) {
      printf DATA "%12.6E ",$period; # transmitter
      printf DATA "%s %8.3f %8.3f ",$site,$lat,$lon; # receiver
      printf DATA "%12.3f %12.3f %12.3f ",$northing{$site},$easting{$site},0.0; # receiver x,y,z
      printf DATA "%s %15.6E %15.6E %15.6E\n",$comp[$i],$dat[2*$i-1],$dat[2*$i],$err[2*$i]; # data
    }
  };
}
close (DATA);

# write a ModEM data file (old format)
#open (DATA, "> DSM1.imp");
#print DATA "Description: Synthetic data set DSM1 from Dublin Institute (2008)\n";
#print DATA "Units: Ohm\n";
#print DATA "Sign convention: +1\n\n";
#print DATA "$periods\n";
#my $period;
#foreach $period (@periods){
#  print DATA "$period\t8\t$sites\n";
#  foreach $site (sort @sites){
#    print DATA "\t$northing{$site}";
#  };
#  print DATA "\n";
#  foreach $site (sort @sites){
#    print DATA "\t$easting{$site}";
#  };
#  print DATA "\n";
#  foreach $site (sort @sites){
#    print DATA "\t0";
#  };
#  print DATA "\n$comment2";
#  foreach $site (sort @sites){
#    #my @d = split(/\s+/,$data{$period}{$site});
#    #my @e = split(/\s+/,$error{$period}{$site});
#    #printf DATA "\t$site %12.3E\n",@d;
#    #printf DATA "\t    %12.3E\n",@e;
#    print DATA "\t$site $data{$period}{$site}\n";
#    print DATA "\t    $error{$period}{$site}\n";
#  };
#}
#close (DATA);
