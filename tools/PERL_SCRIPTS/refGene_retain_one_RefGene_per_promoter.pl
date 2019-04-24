#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

my $refseq = $ARGV[0];
if ($ARGV[0] eq "") {
   $refseq = "/Users/olivier/PROGRAMS/ChIPseeqer-1.0-OLD/DATA/refGene.txt.7June2009";
}

my %ORFPROM = ();
my %NM      = ();

open IN, $refseq;
while (my $l = <IN>) {  
  chomp $l;
  my @a = split /\t/, $l, -1;

  next if ($a[2] =~ /\_/);

  my $orf = $a[12];
  my $tss = undef;
  if ($a[3] eq "+") {
    $tss = $a[4];
  } else {
    $tss = $a[5];
  }

  push @{ $ORFPROM{$orf}{$tss} }, $a[1];

  $NM{ $a[1] } = $l;

}
close IN;



foreach my $g (keys(%ORFPROM)) {
  
  foreach my $tss (keys(%{$ORFPROM{$g}})) {
    my $small = Sets::get_smallest_NM($ORFPROM{$g}{$tss});
    #print "$g\t$tss\t$small\t" . join("\t", @{$ORFPROM{$g}{$tss}}) . "\n";
    print "$NM{$small}\n";
  }

}
