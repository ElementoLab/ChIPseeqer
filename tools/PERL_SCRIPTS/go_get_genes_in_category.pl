BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Sets;
use GroupEnrichment;
use strict;

my $go = GroupEnrichment->new;
$go->setGroups($ARGV[1]);


my $ref = $go->getGeneGroup($ARGV[0]);

Sets::printSet($ref); print "\n";
