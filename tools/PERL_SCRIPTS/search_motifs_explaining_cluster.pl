# for each cluster, find the motifs whose best significant intersection is that cluster
use lib qw(/home/olly/PERL_MODULES);
use ScanACE;
use Table;
use Sets;
use Hypergeom;

my $a_ref_clusterfiles = Sets::getFiles("$ARGV[0]");
my $a_ref_motiffiles   = Sets::getFiles("$ARGV[1]");


my $n = $ARGV[2];

#
#  for each cluster
#

foreach my $r1 (@$a_ref_clusterfiles) {

    print "Cluster $r1\n";

    # get the set of genes
    my $a_ref_set = Sets::readSet($r1);

    my $ta = Table->new;
    $ta->loadFile($r1);

    my $a_ref_set = $ta->getColumn(0);

    #
    #  for each motif, calculate a p-value
    #

    my @a_mots = ();
    
    foreach my $r2 (@$a_ref_motiffiles) {
	
	print " mot $r2\n";
	
	$s = ScanACE->new;
	my $a_ref = $s->_readScanaceResults($r2);
	
	my @a_set = ();
	foreach my $a (@$a_ref) {
	    push @a_set, $a->[0] if (!Sets::in_array($a->[0], @a_set));
	}
	
	#
	#  get the overlap between the two
	#
	my $a_ovl = Sets::getOverlapSet($a_ref_set, \@a_set);
	

	my $s1 = scalar(@$a_ref_set);
	my $s2 = scalar(@a_set);
	my $ov = scalar(@$a_ovl);

	my $pv = Hypergeom::lcumhyper($ov, $s1, $s2, $n);

	my @a_tmp = ($r2, $pv, $ov, $s1, $s2);

	if ($pv < log(0.01)) {
	    push @a_mots, \@a_tmp;
	}
    }

    #
    # order the motifs
    # 
    my @a_mots_sorted = sort { $a->[1] <=> $b->[1] } @a_mots;

    my $i = 0;
    foreach my $r (@a_mots_sorted) {
	print "$r->[0]\t$r->[1]\t$r->[2]\t$r->[3]\t$r->[4]\n";
	$i ++;

	last if ($i == 10);
    }
}

