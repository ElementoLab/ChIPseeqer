#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;



my $nm2orf = "$ENV{CHIPSEEQERDIR}/DATA/refGene.txt.25Nov2009.NM2ORF";

my %ORF2NM = ();
open IN, $nm2orf or die "Cannot open file $nm2orf";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{$ORF2NM{$a[1]}}, $a[0] if (!Sets::in_array($a[0], @{$ORF2NM{$a[1]}}));
}
close IN;




open IN, $ARGV[0] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  print join("\n", @{$ORF2NM{$a[0]}}) . "\n";
}
close IN;

