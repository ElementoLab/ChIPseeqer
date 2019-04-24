use lib qw(/home/olly/PERL_MODULES);

use GO_categories;
use Sets;
use DataFiles;
use GO_func;
use strict;

my $df  = DataFiles->new;

my $m  = GO_func->new;


$m->setSource("GO", "WORM");
$m->setTotalNbORFS($df->get("ELE_BRI")->{"NBGO"});
#$m->setPvalueThreshold(0.001);

$m->setORFset(Sets::readSet($ARGV[0], 1));

my $a_ref_func = $m->getFunctionalContent;
		   
foreach my $r (@$a_ref_func) {
    print sprintf("%3.2e\t%d\t%d\t%d\t%s:%s\n",$r->{PVALUE}, $r->{OVERLAP}, $r->{S1}, $r->{S2}, $r->{NUM}, $r->{TEXT});
}
