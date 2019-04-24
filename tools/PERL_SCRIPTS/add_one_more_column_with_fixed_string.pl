#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";


open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @a, $ARGV[1];
  print join("\t", @a) . "\n";
}
close IN;

