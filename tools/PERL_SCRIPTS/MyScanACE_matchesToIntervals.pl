#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --motifmatches=FILE --ext=INT --showaff=INT\n";
}

my $mot = undef;
my $ext = 0;
my $saf = 0;
my $showall = 0;

GetOptions("motifmatches=s" => \$mot,
	   "showaff=s"      => \$saf,
	   "showall=s"      => \$showall,
           "ext=s"          => \$ext);

open IN, $mot or die "Cannot open $mot\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $c2 = $a[1] + length($a[4]);
  my $c1 = $a[1];
  $c1 = Sets::max(0, $c1-$ext);
  $c2 += $ext;
  print "$a[0]\t$c1\t$c2";
  if ($saf == 1) {
    print "\t$a[3]";
  } elsif ($showall == 1) {
    print "\t" . join("\t", @a[2..$#a]);
  }
  print "\n";
}
close IN;

