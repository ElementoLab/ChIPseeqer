#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;

my $refseq = $ARGV[0];
if ($refseq eq "") {
  $refseq = "/Users/olivier/PROGRAMS/ChIPseeqer-1.0-OLD/DATA/refGene.txt.7June2009";
}

my %ORFPROM = ();

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
  

}
close IN;



foreach my $g (keys(%ORFPROM)) {
  
  foreach my $tss (keys(%{$ORFPROM{$g}})) {
    my $small = Sets::get_smallest_NM($ORFPROM{$g}{$tss});
    print "$g\t$tss\t$small\t" . join("\t", @{$ORFPROM{$g}{$tss}}) . "\n";
  }

}
