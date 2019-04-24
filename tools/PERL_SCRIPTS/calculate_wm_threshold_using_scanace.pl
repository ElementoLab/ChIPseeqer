#
#  using a given WM, run ScanACE, and chose the T that
#   maximizes the overlap with the group of targets ..
#   output threshold and p-value
#
#


BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES /home/elemento/usr/lib/perl5/site_perl/5.6.1/i386-linux";

use ScanACE;
use Motif;
use Hypergeom;


#my $mo = Motif->new;
#$mo->readScanACEMotif($ARGV[0]); # motif

my $sc = ScanACE->new;

$sc->setNbMotifs(5000);
$sc->setMotif($ARGV[0]);
$sc->setFasta($ARGV[1]);
$sc->run;


$sc->setVerbose(1);

my $n     = $ARGV[2];

my $a_ref_set = Sets::readSet($ARGV[3]);

# get ORF -> best score

for (my $nn=1; $nn<=$sc->getNbMotifsInFile($ARGV[0]); $nn++) {

    print ">Motif $nn\n";

    my $a_ref = $sc->getBestSites($nn);
    
    my @SET = ();

    my $sco = 100000000.0;
    my $best_sco = undef;
    my $best_p   = 100000000.0;
    my $best_i   = undef;
    my $best_s1  = undef;
    my $best_s2  = undef;

    foreach my $r (@$a_ref) {
	
	# did we change the threshold ?
	if ($r->[4] < $sco) {
	    
	    # evaluate overlap
	    my $s1 = scalar(@SET);
	    my $s2 = scalar(@$a_ref_set);
	    my $i  = Sets::getOverlapSize(\@SET, $a_ref_set);
	    
	    #Sets::printSet(\@SET);
	    #Sets::printSet($a_ref_set);

	    
	    my $p  = Hypergeom::cumhyper($i, $s1, $s2, $n);
	    
	    #print sprintf("P($i,$s1,$s2,$n)=%3.2e\n", $p);
	    
	    if ($p < $best_p) {
		$best_sco = $sco;
		$best_p   = $p;
		$best_s1   = $s1;
		$best_s2   = $s2;
		$best_i    = $i;
		
	    }

	    $sco = $r->[4];
	} else {
	    
	    # just add to the list
	    push @SET, $r->[0];
	}
	
    }

    my $sp = sprintf("%3.2f", $best_p);
    print "Motif $nn\t$best_i,$best_s1,$best_s2\t$best_sco\t$best_p\n";
    
}


$sc->clean;
