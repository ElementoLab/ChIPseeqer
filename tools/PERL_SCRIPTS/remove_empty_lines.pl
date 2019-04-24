#!/usr/bin/perl

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  if ($l ne "") {
    print "$l\n";
  }
}
close IN;
