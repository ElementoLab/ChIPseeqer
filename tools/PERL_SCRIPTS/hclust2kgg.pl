#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Getopt::Long;
use strict;

my $cdtfile       = undef;
my $gtrfile       = undef;
my $nummaxclust   = 10;

# handling missing arguments
if (@ARGV == 0) {
	die "Usage: perl hlust2kgg --cdt=FILE --gtr=FIRE --clusters=INT\n";
}

GetOptions("cdt=s" => \$cdtfile,
"gtr=s"			=> \$gtrfile,
"clusters=s"    => \$nummaxclust	);

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
    print "GENE\tEXP\n";
    foreach my $k (keys(%SETS)) {
      if (defined($SETS{$k})) {
	foreach my $r (@{$SETS{$k}}) {
	  print "$r\t$p\n";
	}
	$p ++;
      }

    }

    last;
  }

}
close IN;


