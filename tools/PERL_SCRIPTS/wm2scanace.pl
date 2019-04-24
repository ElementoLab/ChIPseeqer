BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Motif;

my $mo = Motif->new;
$mo->readPollardWM($ARGV[0]);


print $mo->generateScanACEMotif;
