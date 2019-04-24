use lib qw(/home/olly/PERL_MODULES);
use Sets;
use GO_func;

$m = GO_func->new;
$m->setSource("GO", "WORM");
$m->setTotalNbORFS(12832);
$m->setPvalueThreshold(0.00001);


my $a_ref = Sets::readSet($ARGV[0], 1);

foreach my $g (@$a_ref) {
    my $a_go_ref = $m->getMIPScategory($g);
    print "$g\n";
    foreach my $f (@$a_go_ref) {
	
	print "\t$f\n";

    }

}

$m->setORFset($a_ref);
$m->setVerbose(0);
my $a_ref_func = $m->getFunctionalContent;
if (scalar(@$a_ref_func) > 0) {
    foreach my $r (@$a_ref_func) {
	my $s_func = lc($r->{TEXT});
	my $s_prob = sprintf("$s_func\t%3.2e", $r->{PVALUE});
	print "$s_prob\n";
    } 
} else {
    print "Nothing\n";
}
