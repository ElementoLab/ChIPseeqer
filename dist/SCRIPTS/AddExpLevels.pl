#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


my $f = "$ENV{HOME}/PROGRAMS/SNPseeqer/TEST/MELNICKRNASEQ/LY1_LY7_CB_NB_RPKM.txt.ORFs";
my %G = ();
open IN, $f or die "Cannot open $f\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $G{$a[5]} = $a[1];
}
close IN;


open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my @b = split /\_/, $l, -1;
  my $e = $G{uc($b[0])};
  if (defined($e)) {
    print "$l\t$e\n"; 
  }
}
close IN;


