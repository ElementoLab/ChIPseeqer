#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

if (@ARGV == 0) {
  die "Args: set orthologs matches motif\n";
}

# read in orthologs
my $h_ref_o = {};
open IN, $ARGV[1] or die "cannot open $ARGV[1]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  #$a[1] =~ s/Pv/PVX_/;

  $h_ref_o->{$a[1]} = $a[0];


}
close IN;


# read in which pf genes have which motifs
open IN, $ARGV[2] or die "cannot open $ARGV[2]\n";
my %M = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $M{$a[0]} = \@a;
}
close IN;


# set
open IN, $ARGV[0] or die "cannot open $ARGV[0]\n";
my $cnto = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  #print "$a[0]\n";
  
  if (defined($h_ref_o->{$a[0]})) {
    
    if (defined($M{$h_ref_o->{$a[0]}})) {

      #print "$a[0]\t" . join("\t", @{$M{$h_ref_o->{$a[0]}}}) . "\n";
    }
    $cnto++;
  }

}
close IN;

print STDERR "$ARGV[3]\t$cnto orthologs\n";


