#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;

my $a_ref1 = [];
open IN, $ARGV[0] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @$a_ref1, "$a[0]\t$a[1]\t$a[2]";
}
close IN;


my $a_ref2 = [];
open IN, $ARGV[1] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @$a_ref2, "$a[0]\t$a[1]\t$a[2]";
}
close IN;



my $a_int = Sets::getOverlapSet($a_ref1, $a_ref2);

Sets::printSet($a_int);
