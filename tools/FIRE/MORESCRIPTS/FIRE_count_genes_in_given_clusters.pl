#
# takes as input list of clusters
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
my $h_ref_p = $ta->getIndexKV(1,0);

# MOTIF = $ARGV[3]

my @GENES = keys(%$h_ref_e);

foreach my $g (@GENES) {
  
  #print "$g: ";
  
  #print $h_ref_p->{ $g };

  if (Sets::in_array($h_ref_e->{ $g }, @a_clusters ) && ($h_ref_p->{$g} eq $ARGV[3])) {
    print "$g\n";
  }
  

}



