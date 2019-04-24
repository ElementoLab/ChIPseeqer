#
# takes as input summary file (with over-rep clusters)
#                signif rep file (with fractions)
#

use lib qw(/home/elemento/PERL_MODULES);

use Fire;
use Sets;
use strict;


# summary
my @a_clusters = split /\,/, $ARGV[0];

# expfile
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref_e   = $ta->getIndexKV(0,1);

# profile
$ta->loadFile($ARGV[2]);
my $a_ref_p = $ta->getArray();

# MOTIF = $ARGV[3]


# 3'UTR lengths
$ta->loadFile($ARGV[4]);
my $h_ref_l = $ta->getIndexKV(0,1);

#print scalar( @{$h_ref_s->{$ARGV[3]}->{CLU}}  ) . " cluster\n";
#print join("-", @{$h_ref_s->{$ARGV[3]}->{CLU}} ) . "\n";

foreach my $g (@$a_ref_p) {


  if (Sets::in_array($h_ref_e->{ $g->[1] }, @a_clusters ) 
      && ($g->[0] eq $ARGV[3])) {
    
    my $p = sprintf("%4.3f", $g->[2] / $h_ref_l->{ $g->[1] });
    
    print "$p\n";
  }


}



