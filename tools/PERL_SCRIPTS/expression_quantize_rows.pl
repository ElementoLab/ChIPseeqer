#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;
print $l;

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $n = shift @a;
  print "$n";
  foreach my $r (@a) {
    if ($r > 1) {
      print "\t2";
    } elsif ($r < -1) {
      print "\t0";
    } else {
      print "\t1";
    }    
  }
  print "\n";
}
close IN;

