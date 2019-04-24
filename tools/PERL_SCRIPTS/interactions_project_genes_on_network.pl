use lib qw(/home/elemento/PERL_MODULES);
use Table;
use Sets;
use strict;

#
# load genes
#
my $h_ref_gen = Sets::getIndex($ARGV[0]);

#
# load interactions
#
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref_int = $ta->getArray();



foreach my $i (@$a_ref_int) {

  if (defined($h_ref_gen->{ $i->[0] }) && defined($h_ref_gen->{ $i->[1] })) {
    print "$i->[0]\t$i->[1]\n";
  }

}

