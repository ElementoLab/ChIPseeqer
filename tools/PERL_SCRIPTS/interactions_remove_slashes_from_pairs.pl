#
#  remove pairs and duplicate the cg matched by the same ps
#

use lib qw(/home/elemento/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  
  my @a1 = split /\//, $r->[0];
  my @a2 = split /\//, $r->[1];
  
  foreach my $c1 (@a1) {
    foreach my $c2 (@a2) {
      print "$c1\t$c2\t$r->[2]\n";
    }
  }
  
}
