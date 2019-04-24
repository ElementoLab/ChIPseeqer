#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

my $f = shift @ARGV;

open IN, $f;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (Sets::in_array($a[1], @ARGV)) {
    $a[1] = join("/", @ARGV);
  } 
  
  print join("\t", @a) . "\n";
  
}
close IN;


