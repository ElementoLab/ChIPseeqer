#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

# read in orthologs
my $h_ref_o = {};
my $h_haso  = {};
open IN, $ARGV[2];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $a[1] =~ s/Pv/PVX_/;

  $h_ref_o->{$a[1]} = $a[0];

  $h_haso->{$a[0]} = $a[1];

  #print "h_ref_o->{$a[1]} = $a[0];\n";

}
close IN;


# open prof 2
my $h_ref_orth = {};
open IN, $ARGV[1];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  if (($a[1] > 0) && (defined($h_ref_o->{$a[0]}))) {
    $h_ref_orth->{  $h_ref_o->{$a[0]}  } = $a[0];
    #print "h_ref_orth->{  $h_ref_o->{$a[0]}  } = 1\n";
  }
  
}
close IN;



# open prof 1

open IN, $ARGV[0];
my $cnttgt = 0;
my $cnttgt_cons = 0;
my $cnttgt_witho = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if ($a[1] > 0) {    
    $cnttgt ++;
    if (defined($h_haso->{$a[0]})) {
      $cnttgt_witho ++;
    }
    if ( defined($h_ref_orth->{ $a[0] } )) {
      #print "$a[0] target and target in ortholog ($h_ref_orth->{$a[0]}).\n";
      $cnttgt_cons++;
    }

  }
}
close IN;
exit if ($cnttgt == 0);
my $frac = sprintf("%4.3f", $cnttgt_cons/$cnttgt_witho);
$ARGV[3] =~ s/\.txt//;
print "$cnttgt\t$cnttgt_witho\t$cnttgt_cons\t$frac\t$ARGV[3]\n";


