#
# INPUT oldrank, oldlist, newlist
#

use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $oldkmer = $a_ref->[ $ARGV[2]-1 ]->[0];


$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
my $rank = 1;
foreach my $r (@$a_ref) {
    if ($r->[0] eq $oldkmer) {
	print "$r->[0]\t$rank\n";
    }

    $rank  ++;
}


