#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


# load genes
my $h_ref_genes = {};
open IN, $ARGV[0] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $h_ref_genes->{$a[0]} = 1;
}
close IN;


#Sets::getIndex($ARGV[0]);

# go thru refgene
my %HE = ();
open IN, $ARGV[1] or die "Cannot open file $ARGV[1]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $n    = $a[1];
  my $c    = $a[2];
  my $fr   = $a[3];
  my $e_st = $a[9];
  my $e_en = $a[10];
  my $e_fr = $a[15];
  my $g    = $a[12];

  next if (!defined($h_ref_genes->{$g}));

  $e_st =~ s/\,$//;
  $e_en =~ s/\,$//;
  $e_fr =~ s/\,$//;

  my @a_e_st = split /\,/, $e_st;
  my @a_e_en = split /\,/, $e_en;
  my @a_e_fr = split /\,/, $e_fr;

  my @exons = ();
  for (my $i=0; $i<@a_e_st; $i++) {
    my @a_tmp = ($a_e_st[$i], $a_e_en[$i]);
    push @exons, \@a_tmp;

  }


  foreach my $e (@exons) {
    my $txte = "$e->[0]\t$e->[1]";
    $HE{$g}{$txte} = 1;
  }


}

close IN;

foreach my $g (keys(%HE)) {
  
  my $sum = 0;
  foreach my $v (keys(%{$HE{$g}})) {
    #print "$v\n";
    my @a = split /\t/, $v;
    my $l = abs($a[1] - $a[0]);
    $sum += $l;
  }
  print "$g\t$sum\n";
}
