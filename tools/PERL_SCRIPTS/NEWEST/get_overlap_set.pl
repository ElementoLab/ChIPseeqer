use lib qw(/home/olly/PERL_MODULES);

use Sets;


my $a_ref1 = Sets::readSet($ARGV[0]);
my $a_ref2 = Sets::readSet($ARGV[1]);


Sets::printSet(Sets::getOverlapSet($a_ref1, $a_ref2));
