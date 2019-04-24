#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Getopt::Long;
use Sets;
use strict;

my $scanacefile = undef;
my $genelist    = undef;
my $n           = 300;
my $namepattern = undef;

if (@ARGV == 0) {
  die "Args: --scanacefile=s --genelist=s --n=s\n";
}

GetOptions("scanacefile=s"    =>  \$scanacefile,
	   "genelist=s"       =>  \$genelist,
	   "n=s"              =>  \$n,
	   "namepattern=s"    =>  \$namepattern);


my $ta = Table->new;

# load genes
my $h_ref_genes = {};
if (defined($genelist)) {
  #print STDERR "Loading genelist ... ";
  $h_ref_genes = Sets::getIndex($genelist);
  #print STDERR "Done.\n";
}




# motif matches
my %MATCHES = ();
$ta->loadFile($scanacefile);
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
  if (!defined($MATCHES{ $r->[0] }) || ($r->[3] > $MATCHES{$r->[0]})) {
    $MATCHES{$r->[0]} = $r->[3];
  }
}

my @genes   = keys(%MATCHES);
my @profile = ();
foreach my $g (@genes) {
  my @a = ($g, 0);
  if (defined($MATCHES{$a[0]}) && ($MATCHES{$a[0]} > 0)) {
    $a[1] = $MATCHES{$a[0]};
  }
  push @profile, \@a; 
}


my @oo = sort { $b->[1] <=> $a->[1] } @profile;

my $t = $oo[$n]->[1];

foreach my $r (@$a_ref) {
  next if (defined($genelist) && !defined($h_ref_genes->{$r->[0]}));
  if ($r->[3] >= $t) {
    if (defined($namepattern)) {
      my $na = undef;
      ($na) = $scanacefile =~ /$namepattern/;
      print "$na\t";
    }
    print Sets::jointab($r);
  }
}


