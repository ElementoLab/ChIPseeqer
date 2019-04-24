#
# making motif target lists
#
# Most of this script was written by Danny Lieber !
#


#!/usr/bin/perl

use lib "$ENV{FIREDIR}/SCRIPTS";

#
# 1. get motif-bin pvalue information as input
# 2. extract gene set for a motif
#


use Sets;
use GroupEnrichment;
use Getopt::Long;
use Table;

use strict;


my $re             = undef;
my $fastafile      = undef;
my $goindex        = undef;
my $gonames        = undef;
my $maxcatsize     = 500;
my $pmin           = 0.01;

my $threshold      = 0.05; # pvalue threshold for bin significance 
my $bonf_threshold =  undef; # bonferroni-corrected threshold  

my $outgo          = undef;
my $outfile        = undef;

my $profiles       = undef;
my $expfile        = undef;
my $matrixfile     = undef;
my $clusterfile    = undef;
my $usemodule      = undef;
my $nbcols         = undef;
my $nbrows         = undef;
my @temparray      = undef;

if (@ARGV == 0) {
  die "Usage: perl get_motif_targets.pl --matrixfile=FILE --expfile=FILE --profiles=FILE --outfile=FILE\n";
}

GetOptions ('matrixfile=s'     => \$matrixfile,
	    'profiles=s'       => \$profiles,
	    'expfile=s'        => \$expfile,
            'outfile=s'        => \$outfile

	  );


die "Please specify outfile .. \n" if (!defined($outfile));

use Table;

my $ta = Table->new;


#
# read in FULL over-rep p-value matrix (must be the full one)
#
my %CLU = ();
my @MOTIFS = ();
$ta->loadFile($matrixfile);
my $a_ref_sum = $ta->getArray();
shift @$a_ref_sum;


# get dimensions
$nbcols = $ta->getNbColumns;
$nbrows = $ta->getNbRows;

# get threshold for over-rep
$bonf_threshold = $threshold / $nbcols;
my $logbonf = - Sets::log10($bonf_threshold);


# get set of clusters where gene is overrep
foreach my $r (@$a_ref_sum) {
  push @MOTIFS, $r->[0];
  for (my $i=1; $i<@$r; $i++) {
    if ($r->[$i] > $logbonf) {
      $CLU{$r->[0]}{$i-1} = 1; # bin i-1 is at column i
    }
  }
}

#
#  read in expfile
#
$ta->loadFile($expfile);
my $h_ref_exp = $ta->getIndexKV(0,1);
my $a_ref_exp = $ta->getColumn(0);

#
# load profiles
#
my %CLU_GENES = ();
my %CLU_GENES_ALL = ();
$ta->loadFile($profiles);
my $a_ref_pro = $ta->getArray();

# loop over motif occurrences
foreach my $r (@$a_ref_pro) {

  # get the expression cluster 
  my $c = $h_ref_exp->{ $r->[1] };
  next if (!defined($c));

  # get the motif 
  my $m = $r->[0];

  # add to the ALL list
  push @{ $CLU_GENES_ALL{ $m }}, $r->[1];

  # add only if belongs to overexp cluster
  if (defined($CLU{$m}{$c})) {
    push @{ $CLU_GENES{ $m }}, $r->[1];
  }
}


open OUT, ">$outfile" or die "Cannot open $outfile.\n";
open OUT_ALL, ">$outfile.all" or die "Cannot open $outfile.all.\n";

foreach my $k (@MOTIFS) {

  my $s  = Sets::removeDuplicates($CLU_GENES{$k});
  
  print OUT $k, "\t";
  print OUT join "\t", @{$s};
  print OUT "\n";
  
  $s  = Sets::removeDuplicates($CLU_GENES_ALL{$k});
  
  print OUT_ALL $k, "\t";
  print OUT_ALL join "\t", @{$s};
  print OUT_ALL "\n";
  
  
  
}

close OUT;
