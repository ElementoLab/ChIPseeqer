use lib qw(/home/olly/PERL_MODULES);

use GO_categories;
use Sets;

my $goc = GO_categories->new;



#$goc->setVerbose(1);
$goc->setID("root", "", "localhost");
$goc->setSpecies("HUMAN");

$goc->setTotalNbORFS(8706);

if ($ARGV[1]) {
    $goc->setEnsemblORFset(Sets::readSet($ARGV[0], 1));
} else {
    $goc->setORFset(Sets::readSet($ARGV[0], 1));
}
#$goc->setBonferroni(0);
$goc->setPvalueThreshold(1.0);


my $a_ref = $goc->getFunctionalContent;

foreach my $r (@$a_ref) {
    print sprintf("%3.2e\t%d\t%d\t%d\t%s:%s\n",$r->{PVALUE}, $r->{OVERLAP}, $r->{S1}, $r->{S2}, $r->{NUM}, $r->{TEXT});
    
}
