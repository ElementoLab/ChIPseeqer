#!/usr/bin/perl

open IN, $ARGV[0];

my $l = <IN>;

while (my $l = <IN>) {
  print $l;
}
close IN;

