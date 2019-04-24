use lib qw(/home/olly/PERL_MODULES);
use Sets;
use COG;

my $co = COG->new;

$co->setSpecies($ARGV[1]);
$co->setTotalNbORFS($ARGV[2]);
$co->setBonferroni(1);
open IN, $ARGV[0];
my $pheno_single = 0;
while (my $l = <IN>) {
    chomp $l;

    if ($l =~ /^PHENOTYPE (.+?)\t(.+)$/) {
        my $t = $2;
	$p = $1;

	print "$lastp\n"; 

        if ($t =~ /GENES/) {
            $pheno_single = 1;
	    
	    
	    
	    if (scalar(@bag) > 0) {
		#
		#  
		#
		$co->setORFset(\@bag);
		my $a_ref = $co->getFunctionalContent();
		foreach my $c (@$a_ref) {
		    print "\t$c->{TEXT}\t" . sprintf("%3.2e\n", $c->{PVALUE});
		}
		@bag = ();
	    }

	    

        } else {
            exit;
        }

	$lastp = $p;
    } else {

        if ($pheno_single == 1) {
            
            if ($l =~ /^>([^>].+?)\t([\d\.]+?)\t(.+?)\t/) {
		
		my $s = $2; next if ($2 < 0.1);
                #print "$1\t$3\n"; 
		#next if (scalar(@bag) == 50);
		push @bag, $3;

            }

        }


    }


}



