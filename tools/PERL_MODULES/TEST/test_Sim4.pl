BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sim4;

my $si = Sim4->new;

$si->run($ARGV[0], $ARGV[1]);


