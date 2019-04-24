use lib qw(/home/elemento/PERL_MODULES);
srand();

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

$ta->randomizeColumn($ARGV[0]);

$ta->print;
