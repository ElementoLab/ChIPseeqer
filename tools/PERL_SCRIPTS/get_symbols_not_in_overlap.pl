use lib qw(/Users/olivier/PERL_MODULES);
use Sets;

my $a_ref1 = Sets::readSet($ARGV[0]);
my $a_ref2 = Sets::readSet($ARGV[1]);

my $a_not1 = Sets::getNotIntersection($a_ref1, $a_ref2, 0);
my $a_not2 = Sets::getNotIntersection($a_ref1, $a_ref2, 1);


print "In 1 but not in 2:\n";
Sets::printSet($a_not1);
print "\n";
print "In 2 but not in 1:\n";
Sets::printSet($a_not2);

