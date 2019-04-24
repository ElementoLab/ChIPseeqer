use lib qw(/home/olly/PERL_MODULES);

use GO_categories;
use Sets;
use DataFiles;

my $goc = GO_categories->new;
my $df  = DataFiles->new;


$goc->setID("root", "", "localhost");
$goc->setSpecies($ARGV[1]);

$goc->setTotalNbORFS($df->get("ELE_BRI")->{NBGO});

#if ($ARGV[1]) {
#    $goc->setEnsemblORFset(Sets::readSet($ARGV[0], 1));
#} else {
$goc->setORFset(Sets::readSet($ARGV[0], 1));
#}

#$goc->setPvalueThreshold(0.0001);


my $a_ref = $goc->getFunctionalContent;

foreach my $r (@$a_ref) {
    print sprintf("%3.2e\t%d\t%d\t%d\t%s:%s\n",$r->{PVALUE}, $r->{OVERLAP}, $r->{S1}, $r->{S2}, $r->{NUM}, $r->{TEXT});
}

