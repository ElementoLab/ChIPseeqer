BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
shift @$a_ref;
my $t = $ARGV[1];

print "GENE\tEXP\n";
foreach my $r (@$a_ref) {

  if (($r->[1] < $t) && ($r->[2] < 0)) {
    print "$r->[0]\t1\n";
  } else {
    print "$r->[0]\t0\n";
  }

}

