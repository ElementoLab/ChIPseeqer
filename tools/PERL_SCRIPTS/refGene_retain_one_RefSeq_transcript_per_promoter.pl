#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;

my $refseq = $ARGV[0];
#if ($ARGV[0] eq "") {
#   $refseq = "/Users/olivier/PROGRAMS/ChIPseeqer-1.0-OLD/DATA/refGene.txt.7June2009";
#}

my %ORFPROM = ();

open IN, $refseq;
my %NM = ();
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
  $NM{$a[1]} = \@a;

}
close IN;



foreach my $g (keys(%ORFPROM)) {
  
  foreach my $tss (keys(%{$ORFPROM{$g}})) {
    my $small = Sets::get_smallest_NM($ORFPROM{$g}{$tss});
    #print "$g\t$tss\t$small\t" . join("\t", @{$ORFPROM{$g}{$tss}}) . "\n";
    if ($ARGV[1] eq "") {
      print join("\t", @{$NM{$small}}) . "\n";
    } else {
      my @a =  @{$NM{$small}};
      my $st = undef;
      if ($a[3] eq '+') {
	$st = 1;
      } else {
	$st = -1;
      }
      print "$a[1]\t$a[2]\t$a[4]\t$a[5]\t$st\n";
    }
  }

}
