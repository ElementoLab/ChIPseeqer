#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

my %H = ();
open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $H{$a[0]} = $a[1];
}
close IN;


open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (defined($H{$a[0]}) && ($H{$a[0]} == 1)) {
    print "$a[0]\t$a[1]\n";
  }

}
close IN;


