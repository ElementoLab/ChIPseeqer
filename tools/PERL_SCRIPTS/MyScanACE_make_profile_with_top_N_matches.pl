#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Getopt::Long;
use Sets;
use strict;

my $scanacefile = undef;
my $genelist    = undef;
my $n           = 250;
my $outdir      = undef;

if (@ARGV == 0) {
  die "Args: --scanacefile=s --genelist=s [ or --fastafile ] --n=s\n";
}

GetOptions("scanacefile=s"    =>  \$scanacefile,
	   "genelist=s"       =>  \$genelist,
	   "n=s"              =>  \$n,
	   "outdir=s"         =>  \$outdir);


my $ta = Table->new;

# load genes
$ta->loadFile($genelist);
my $a_ref_genes = $ta->getColumn(0);

# motif matches
my %MATCHES = ();
$ta->loadFile($scanacefile);
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {

  if (!defined($MATCHES{ $r->[0] }) || ($r->[3] > $MATCHES{$r->[0]})) {
    $MATCHES{$r->[0]} = $r->[3];
  }
  
}

my @profile = ();
foreach my $g (@$a_ref_genes) {
  my @a = ($g, 0);
  if (defined($MATCHES{$a[0]}) && ($MATCHES{$a[0]} > 0)) {
    $a[1] = $MATCHES{$a[0]};
  }


  push @profile, \@a; 
}


my @oo = sort { $b->[1] <=> $a->[1] } @profile;

my $i = 0;
print "GENE\tEXP\n";
foreach my $r (@oo) {
  if (($i <= $n) && ($r->[1] > 0)) {
    print "$r->[0]\t1\n";
  } else {
    print "$r->[0]\t0\n";
  }
  $i++;
}
