#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


open IN, $ARGV[0] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $l = abs($a[2] - $a[1]);

  print "$l\n";

}
close IN;

