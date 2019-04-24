use lib qw(/home/olly/PERL_MODULES);

use GenRegexp;

my $ge = GenRegexp->new;
$ge->calc($ARGV[0], $ARGV[1]);

print $ge->getNbMatches(); print "\n";
