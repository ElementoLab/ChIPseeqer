# load .summary (where motifs are overrep)
#      .expfile (which cluster does a gene belong to)
#      .profiles (whether a gene has a motif)

#
# takes as input summary file (with over-rep clusters)
#                signif rep file (with fractions)
#


use lib "$ENV{FIREDIR}/SCRIPTS";

use Fire;
use Sets;
use strict;

my $expfile = $ARGV[0];

# summary
my $a_ref_s = Fire::loadFireMotifSummaryArray("$expfile.summary");

# expfile
my $ta = Table->new;
if (-e "$expfile.quantized") {
  $ta->loadFile("$expfile.quantized");
} else {
  $ta->loadFile($expfile);
}
my $h_ref_e   = $ta->getIndexKV(0,1);

# profiles
$ta->loadFile("$expfile.profiles");
my $a_ref_p = $ta->getArray(); #getIndexKV(1,0);
my $h_ref_p = {};
foreach my $r (@$a_ref_p) {
  push @{ $h_ref_p->{$r->[1]} }, $r->[0];
}

# MOTIF = $ARGV[3]

my @GENES = keys(%$h_ref_e);

# iterate other motifs

open OUT, ">$expfile.targets" or die "Cannot open outfile\n";

foreach my $m (@$a_ref_s) {

  print OUT "$m->{MOTIF}";
  foreach my $g (@GENES) {
    if (Sets::in_array($h_ref_e->{ $g }, @{ $m->{CLU} } ) && Sets::in_array($m->{MOTIF}, @{$h_ref_p->{$g}})) { 
      print OUT "\t$g";
    }
  }
  print OUT "\n";
  
}
close OUT;


