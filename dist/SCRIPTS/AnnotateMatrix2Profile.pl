#!/usr/bin/perl
use strict;

open IN, $ARGV[0];
my $l = <IN>; chomp $l;
my @a = split /\t/, $l;

print "$a[0]\tALL\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $istgt = 0;
  foreach my $p (@a) {
    if ($p > 0) {
      $istgt = 1;
      last;
    }
  }
  print "$a[0]\t$istgt\n";

}
close IN;

