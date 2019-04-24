use lib qw(/home/elemento/PERL_MODULES);
use GeneWise;

my $ge = GeneWise->new;

$ge->_get_result($ARGV[0]);
