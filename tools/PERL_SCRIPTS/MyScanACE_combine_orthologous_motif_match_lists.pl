use strict;

my $orth = "/Users/olivier/PROGRAMS/PLASMODIUM/GENOMES/proteins/pf_pv_orthologs.txt";

open IN, $orth;
my $h_orth = 
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
}
close IN;

my %MATCHES = ();
open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my @a_tmp = ($a
}
close IN;



