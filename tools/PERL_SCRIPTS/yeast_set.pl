#
#  do a complete analysis of a set of yeast ORFs 
#
use lib qw(/home/olly/PERL_MODULES);
use GO_func;
use Sets;
use Yeast;

#  infos about the set to analyze
my $a_ref = Sets::readSet($ARGV[0]);
my $nset  = scalar(@$a_ref);


# infos about the total number of genes
my $ntotal = 6000;


#
#  MIPS categories
#
my $m = GO_func->new;
$m->setSource("MIPS", "YEAST");
$m->setORFset($a_ref);

$m->setTotalNbORFS($ntotal);

$m->setPvalueThreshold(0.0001);

print "-" x 100; print "\n";
print "MIPS (p<0.001) \n";
print "-" x 100; print "\n";

foreach my $r (@{$m->getFunctionalContent}) {
    print "$r->{NUM}\t" . sprintf("%e", $r->{PVALUE}) . "\t" . sprintf("%3d", $r->{OVERLAP}) . " / $r->{TOTAL}\t" . lc($r->{TEXT});   
    print "\n";
    
    #last;
}

print "\n";


#
#  ChIP enrichment
#
print "-" x 100; print "\n";
print "ChIP (Young data, p<0.001) \n";
print "-" x 100; print "\n";

my $y = Yeast->new;
my $a_ref_chip = $y->getCHIPEnrichment($a_ref, $ntotal, 0.0001);

my $cnt = 0;
foreach my $r (@$a_ref_chip) {
    if ($r->{P} < 0.0001) {
	print "$r->{REGULATOR}\t" . sprintf("%e", $r->{P}) . "\n";
	$cnt ++;
    }
}

print "Nothing ...\n" if ($cnt == 0);
print "\n";


#


#
#  motif content
#

print "-" x 100; print "\n";
print "Known motifs\n";
print "-" x 100; print "\n";
my $a_mot = $y->getMotifEnrichment($a_ref, $ntotal, 0.0001);
foreach my $m (@$a_mot) {
    print sprintf("$m->[0]\t%e\n", $m->[1]);
}


exit;

#  RAP1 enrichment
#
my $a_ref_brown_RAP1 = Sets::readSet("/home/olly/DATA/CHIP/BROWN/RAP1.txt");
my $a_ovl_RAP1       = Sets::getOverlapSet($a_ref_brown_RAP1, $a_ref);
my $s1               = scalar(@$a_ref_brown_RAP1);
my $s2               = scalar(@$a_ref);
my $i                = scalar(@$a_ovl_RAP1);
my $p                = Hypergeom::cumhyper($i, $s1, $s2, $ntotal);
print "-" x 100; print "\n";
print "ChIP RAP1 (P. Brown) \n";
print "-" x 100; print "\n";
if ($p < 0.01) {
    print sprintf("%e\n", $p);
} else {
    print "Nothing ...\n";
}
print "\n";



#
#  evolutionary divergence with bayanus
#

my $a_ref_div = $y->getDivergences;
my $a_ref_dis = Sets::SQLRefToSet($a_ref_div, "DISTANCE_JC_BAYANUS");
# get the median distance
my $med = Sets::median($a_ref_dis);
# get the number of 
my $h_ref = Sets::SQLRefToSimpleHash($a_ref_div, "ORF", "DISTANCE_JC_BAYANUS");
my $less = 0; my $more = 0;
foreach my $r (@$a_ref) {
    
    if (defined($h_ref->{$r}) && $h_ref->{$r} < $med) {
	$less ++;
    }

    if (defined($h_ref->{$r}) && $h_ref->{$r} > $med) {
	$more ++;
    }
    
}
my $total = $less + $more;
my $pb = Hypergeom::cumbino($more, $total, 0.5);

print "-" x 100; print "\n";
print "Divergence with Bayanus \n";
print "-" x 100; print "\n";
print "Overall median is $med, $more/$total > median in this set, p=" . sprintf("%e", $pb) . "\n"; 

print "\n";


#
#  mRNA level
#
my $a_ref_div = $y->getYeastFeature("HOSLTEGE_MRNA_LEVEL");
my $a_ref_dis = Sets::SQLRefToSet($a_ref_div, "HOSLTEGE_MRNA_LEVEL");
# get the median distance
my $med = Sets::median($a_ref_dis);
# get the number of 
my $h_ref = Sets::SQLRefToSimpleHash($a_ref_div, "ORF", "HOSLTEGE_MRNA_LEVEL");
my $less = 0; my $more = 0;
foreach my $r (@$a_ref) {
    
    if (defined($h_ref->{$r}) && $h_ref->{$r} < $med) {
	$less ++;
    }

    if (defined($h_ref->{$r}) && $h_ref->{$r} > $med) {
	$more ++;
    }
    
}
my $total = $less + $more;
my $pb = Hypergeom::cumbino($more, $total, 0.5);

print "-" x 100; print "\n";
print "mRNA level\n";
print "-" x 100; print "\n";
print "Overall median is $med, $more/$total > median in this set, p=" . sprintf("%e", $pb) . "\n"; 

print "\n";

#
#  mRNA half-life
#


#
#  mRNA transcrition rate
#


#
#  protein levels
#

my $a_ref_div = $y->getYeastFeature("PROTEIN_ABUNDANCE");
my $a_ref_dis = Sets::SQLRefToSet($a_ref_div, "PROTEIN_ABUNDANCE");
# get the median distance
my $med = Sets::median($a_ref_dis);
# get the number of 
my $h_ref = Sets::SQLRefToSimpleHash($a_ref_div, "ORF", "PROTEIN_ABUNDANCE");
my $less = 0; my $more = 0;
foreach my $r (@$a_ref) {
    
    if (defined($h_ref->{$r}) && $h_ref->{$r} < $med) {
	$less ++;
    }

    if (defined($h_ref->{$r}) && $h_ref->{$r} > $med) {
	$more ++;
    }
    
}
my $total = $less + $more;
my $pb = Hypergeom::cumbino($more, $total, 0.5);

print "-" x 100; print "\n";
print "Protein abundance\n";
print "-" x 100; print "\n";
print "Overall median is $med, $more/$total > median in this set, p=" . sprintf("%e", $pb) . "\n"; 

print "\n";


goto TOTO;
#
#  co-expression Gasch + Spellman
#

# create a ORF => nb for faster data access
my %hash = ();
my $i = 0;
foreach my $r (@$a_ref) {
    $hash{$r} = $i;
    $i++;
}

# read the entire correlation file and get an array of corrs for the set, and a global array of corrs
open IN, "/home/olly/PROGRAMS/COEXPRESSION/GASP_CORRELATIONS.txt";
my @a_loc  = ();
my @a_glo  = ();
my @a = ();
while (my $l = <IN>) {
    chomp $l;
    
    @a = split /\t/, $l;
    if (defined($hash{$a[0]}) && defined($hash{$a[1]})) {
	push @a_loc, $a[2];
    }
    push @a_glo, $a[2];
}
close IN;


# calculate the median
print "Calc median ..";
my $med = Sets::median(\@a_glo);
print "done\n";

# get the number of correlation values above and below the median

my $less = 0; my $more = 0;
foreach my $r (@a_loc) {
    
    if ($r < $med) {
	$less ++;
    }

    if ($r > $med) {
	$more ++;
    }
    
}
my $total = $less + $more;
my $pb = Hypergeom::cumbino($more, $total, 0.5);

print "-" x 100; print "\n";
print "Expression correlation\n";
print "-" x 100; print "\n";
print "Overall median is $med, $more/$total > median in this set, p=" . sprintf("%e", $pb) . "\n"; 

print "\n";

TOTO:
#
#  co-expression Spellman
#


#
#  co-expression Rosetta
#



#
#  essentiality
#
my $a_ess = $y->getEssentialGenes();
my $p     = scalar(@$a_ess) / $ntotal;
my $a_ess_in_set = Sets::getOverlapSet($a_ess, $a_ref);
my $k     = scalar(@$a_ess_in_set);
my $N     = scalar(@$a_ref);

print "-" x 100; print "\n";
print "Essentiality\n";
print "-" x 100; print "\n";
my $pb = Hypergeom::cumbino($k, $N, $p);
print "$k out of $N genes are essential in this set, p=" . sprintf("%e", $pb). "\n";
print "\n";
print "\n";
