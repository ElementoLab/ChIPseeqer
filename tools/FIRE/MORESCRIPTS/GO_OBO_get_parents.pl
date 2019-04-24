BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use GO;

my $go = GO->new;

$go->load_OBO_Ontology($ARGV[0]);
$go->printNodes();
