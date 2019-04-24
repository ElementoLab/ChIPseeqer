#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --f1=FILE --f2=FILE\n";
}


my $motifmatches = undef;
my $peakfile     = undef;

GetOptions("motifmatches=s" => \$motifmatches,
           "peakfile=s"     => \$peakfile);

my $txt = `CompareIntervals -peakfile1 $peakfile -peakfile2 $motifmatches -show_ov_int 1`;

my @a_txt = split /\n/, $txt;

foreach my $l (@a_txt) {
  my @a = split /\t/, $l;
  my $c = shift @a;
  my $i = shift @a;
  my $j = shift @a;
  my $n = shift @a;
  my $l = $j - $i + 1;
  foreach my $t (@a) {
    my @b = split /\-/, $t;
    my $d = $b[2] - $b[1];
    my $p = int(rand($l - $d));  
    $p = $i + $p;
    my $pd = $p+$d;
    print "$c\t$p\t$pd\n";    
  }
}

