#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

open IN, $ARGV[0];

my @M = ();

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (($a[0] =~ /\!Sample_title/) || ($a[0] =~ /\!Sample_characteristics_ch1/)) {
    #print join("\t", @a) . "\n";
    push @M, \@a;
  }
  

}
close IN;


my $a_ref_t = Sets::transpose(\@M);
shift @$a_ref_t;
foreach my $t (@$a_ref_t) {
  #print join("\t", @$t) . "\n";
  $t->[0] =~ s/\"//g;
  $t->[8] =~ s/\"//g;
  $t->[9] =~ s/\"//g;

  # CHOP
  #$t->[10] =~ s/\"//g;
  #if (
  print "$t->[0]\t$t->[8]\t$t->[9]\n";
}
