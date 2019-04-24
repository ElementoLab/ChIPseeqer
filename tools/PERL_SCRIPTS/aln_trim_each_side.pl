
use lib "$ENV{HOME}/PERL_MODULES";

use ClustalW;
use strict;

my $cl = ClustalW->new;
#$cl->setNumCharName(length("031008-E4       "));
$cl->setFile($ARGV[0]);

$cl->removeTrailingSequences($ARGV[1], 0);
$cl->removeTrailingSequences($ARGV[2], 1);

#print $cl->getClustalWformat;

print $cl->getFastaFormat();

