#
# takes as input summary file (with over-rep clusters)
#                signif rep file (with fractions)
#
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

#use lib qw(/home/elemento/PROGRAMS/MIMOTIFS);

use Fire;
use Sets;
use strict;


# $ARGV[0] = motif

# expfile
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref_e   = $ta->getIndexKV(0,1);

# profiles
$ta->loadFile($ARGV[2]);
my $a_ref_p = $ta->getArray();




my %H = ();

foreach my $r (@$a_ref_p) {
  
  if (defined($h_ref_e->{$r->[1]})) {
    push @{ $H{ $h_ref_e->{$r->[1]} } }, $r;
  }
}


foreach my $k (sort { $a <=> $b } keys(%H)) {
  

  foreach my $r ( @{ $H{ $k } }) {
    print "Cluster$k\t";

    print join("\t", @$r); print "\n";
  }
  
}


