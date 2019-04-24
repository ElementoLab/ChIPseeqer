#!/usr/bin/perl

use strict;

open IN, $ARGV[0];
my $l = <IN>;
my @COUNT = ();
while (my $l = <IN>) {

  chomp $l;
  my @a = split /\t/, $l, -1;

  $COUNT[0] ++;
  for (my $i=1; $i<@a; $i++) {
    if ($a[$i] > 0) {
      $COUNT[$i] ++;
    }
  }

}
close IN;


open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my @newa = ();
  for (my $i=0; $i<@a; $i++) {
    if ($COUNT[$i] > $ARGV[1]) {
      push @newa, $a[$i];
    }
  }
  print join("\t", @newa);
  print "\n";
  
}
close IN;

