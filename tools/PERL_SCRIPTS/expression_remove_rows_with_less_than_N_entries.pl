#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $numreal = 0;
  foreach my $e (@a) {
    if ($e ne "") {
      $numreal++;
    }
  }
  if ($numreal >= $ARGV[1]) {
    print join("\t", @a) . "\n";
  }
  
}
close IN;

