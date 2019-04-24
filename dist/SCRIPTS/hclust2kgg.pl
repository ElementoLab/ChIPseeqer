#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}}";
use Getopt::Long;
use strict;

my $cdtfile       = undef;
my $gtrfile       = undef;
my $nummaxclust   = 10;
my $refgene       = "$ENV{CHIPSEEQERDIR}/DATA/refGene.txt.07Jun2010";
my $usebkg        = 1;

# handling missing arguments
if (@ARGV == 0) {
  die "Usage: perl hlust2kgg --cdt=FILE --gtr=FILE --clusters=INT\ --refgene=FILE n";
}

GetOptions(
	   "cdt=s"			=> \$cdtfile,
	   "gtr=s"			=> \$gtrfile,
	   "clusters=s"		=> \$nummaxclust,
	   "usebkg=s"		=> \$usebkg,
	   "refgene=s"		=> \$refgene,
	  );

# read refGene first to get the total number of genes
open IN, $refgene or die "Cannot open $refgene\n";
my %REFGENE = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $REFGENE{$a[1]} = 1;
}
close IN;



# how many genes do we have, read CDT
open IN, $cdtfile;
my $l = <IN>;
$l = <IN>;

my $numgenes = 0;
my %SETS = ();

while (my $l = <IN>) {  
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{$SETS{$a[0]}}, $a[1];
  $numgenes ++;

}
close IN;

# now agglomerate gene sets
open IN, $gtrfile;
my $numclust = $numgenes;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  push @{$SETS{$a[0]}}, @{$SETS{$a[1]}};
  push @{$SETS{$a[0]}}, @{$SETS{$a[2]}};

  $SETS{$a[1]} = undef;
  $SETS{$a[2]} = undef;

  $numclust --;
	
  if ($numclust == $nummaxclust) {
    my $p = 0;
    if ($usebkg == 1) {
      $p = 1;
    }	
    print "GENE\tEXP\n";
    foreach my $k (keys(%SETS)) {
      if (defined($SETS{$k})) {
	foreach my $r (@{$SETS{$k}}) {
	  print "$r\t$p\n";
	  if (defined($REFGENE{$r}) && ($REFGENE{$r} == 1)) {
	    $REFGENE{$r} = 0;
	  }
	}
	$p ++;
      }

    }

    last;
  }

}
close IN;

if ($usebkg == 1) {
  foreach my $g (keys(%REFGENE)) {
    if ($REFGENE{$g} == 1) {
      print "$g\t0\n";
    }    
  }
}

