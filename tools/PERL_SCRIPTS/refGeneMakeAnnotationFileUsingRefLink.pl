#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

if (@ARGV == 0) {
  die "Args: refgene reflink\n";
}

# load refgene 
my %T = ();
open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $T{$a[1]} = 1;
}
close IN;


# load reflink 
open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if (defined($T{$a[2]})) {
    print "$a[2]\t$a[0]\t$a[1]\n";
  }
}
close IN;

