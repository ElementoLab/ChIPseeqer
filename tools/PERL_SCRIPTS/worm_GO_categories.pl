use lib qw(/home/olly/PERL_MODULES);
use Sets;
use GO_categories;
use Table;
use strict;
use DataFiles;

my $df      = DataFiles->new;


my $p_global = (defined($ARGV[1])?$ARGV[1]:0.00001);
#my $p_global = 1.0; #0.00001;

my $ntotal   = 12832;
#my $ntotal   = 17468;
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_dup = $ta->getColumn(0);
my $a_ref = Sets::removeDuplicates($a_ref_dup);
my $n  = scalar(@$a_ref);



my $goc = GO_categories->new;
$goc->setVerbose(0);
$goc->setID($df->get("USER"), $df->get("PASS"), $df->get("HOST"));
$goc->setSpecies("WORMS");
#$goc->setOriginal(1);
#$goc->setMaxCategory(1000);
$goc->setMinCategory(5);
$goc->setTotalNbORFS($ntotal);

my $h_ref_cat = $goc->getGOCategoriesNames();
$goc->setPvalueThreshold($p_global);
if ($ARGV[1]) {
foreach my $g (@$a_ref) {

    my $a_ref_cat = $goc->getGOCategories($g);
    print "$g\t";
    foreach my $c (@$a_ref_cat) {
	#print "\t"; 
	print $h_ref_cat->{$c}->{NAME}; 
	print "\t";
    }
   print "\n";
}
}
print "\n--> Functional enrichment (GO) for the $n genes :\n";


$goc->setORFset($a_ref);
$goc->setOriginal(0);

print "1. Biological processes:\n";
$goc->setAnnotationType("P");
my $a_ref_func = $goc->getFunctionalContent;
if (scalar(@$a_ref_func) > 0) {
    
    foreach my $c (@$a_ref_func) {
	my $s_func = lc($c->{TEXT});
	my $s_prob = sprintf("p=%3.2e", $c->{PVALUE});
	print "$s_prob\t$s_func\t($c->{OVERLAP} in this set, $c->{EXP} expected)\n";
    }
    
} else {
    print "Nothing significant\n";
    
}



print "\n2. Molecular function:\n";
$goc->setAnnotationType("F");
my $a_ref_func = $goc->getFunctionalContent;
if (scalar(@$a_ref_func) > 0) {
    
    foreach my $c (@$a_ref_func) {
	my $s_func = lc($c->{TEXT});
	my $s_prob = sprintf("p=%3.2e", $c->{PVALUE});
	print "$s_prob\t$s_func\t($c->{OVERLAP} in this set, $c->{EXP} expected)\n";
    }
    
} else {
    print "Nothing significant\n";
    
}

print "\n3. Cellular localization:\n";
$goc->setAnnotationType("C");
my $a_ref_func = $goc->getFunctionalContent;
if (scalar(@$a_ref_func) > 0) {
    
    foreach my $c (@$a_ref_func) {
	my $s_func = lc($c->{TEXT});
	my $s_prob = sprintf("p=%3.2e", $c->{PVALUE});
	print "$s_prob\t$s_func\t($c->{OVERLAP} in this set, $c->{EXP} expected, $c->{S1} total in GO)\n";
    }
    
} else {
    print "Nothing significant\n";
    
}
