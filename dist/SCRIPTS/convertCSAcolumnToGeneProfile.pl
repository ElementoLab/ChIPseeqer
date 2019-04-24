#!/usr/bin/perl
use lib "$ENV{CHIPSEEQERDIR}";
use Sets;
use strict;

my %H = ();
open IN, "$ENV{CHIPSEEQERDIR}/DATA/hg18/refGene.txt.07Jun2010.NM2ORF" or die "Cannot open $_\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $H{$a[0]} = $a[1];
}
close IN;


open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;
print $l;

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[1] > 1) {
    $a[1] = 1;
  }
  if (defined($H{$a[0]})) {
    print "$H{$a[0]}\t$a[1]\n";
  }
}
close IN;

