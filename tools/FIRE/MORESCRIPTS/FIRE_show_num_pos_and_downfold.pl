#
# takes as input summary file (with over-rep clusters)
#                signif rep file (with fractions)
#

use lib qw(/home/elemento/PROGRAMS/MIMOTIFS);

use Fire;
use Sets;
use strict;


# summary
my $h_ref_s = Fire::loadFireMotifSummary($ARGV[0]);

# expfile
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref_e   = $ta->getIndexKV(0,1);

# profile
$ta->loadFile($ARGV[2]);
my $a_ref_p = $ta->getArray();

# MOTIF = $ARGV[3]

# matrix
$ta->loadFile($ARGV[4]);
my $a_ref_m = $ta->getArray(0);

my %count = ();
foreach my $g (@$a_ref_p) {

  if (Sets::in_array($h_ref_e->{ $g->[1] }, @{$h_ref_s->{$ARGV[3]}->{CLU}} ) 
      && ($g->[0] eq $ARGV[3])) {
    $count{ $g->[1] } ++;
    
  }
}

foreach my $r (@$a_ref_m) {
  if (defined($count{$r->[0]})) {
    print "$count{$r->[0]}\t";
    print join("\t", @$r);
    print "\n";
  }
}



