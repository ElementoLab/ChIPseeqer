use lib qw(/home/olly/PERL_MODULES);
use Sets;

use Fasta;
use strict;
use Hypergeom;


#
#  read in the 7-mers
#
my $a_ref_kmers = Sets::readSet($ARGV[0]);
my @H = ();

my $k = undef;
foreach my $r (@$a_ref_kmers) {
    $H[0]->{ $r } = 0;
    $H[1]->{ $r } = 0;
    $H[2]->{ $r } = 0;
    $k = length($r);
}



#
#  traverse the non-aligned, than the aligned sequence
#
for (my $j=0; $j<2; $j++) {

    my $fa = Fasta->new;
    $fa->setFile($ARGV[$j+1]);
    
    while (my $a_ref = $fa->nextSeq()) {
	my ($n, $s) = @$a_ref; 

	my $l = length($s);
	for (my $i=0; $i<$l-$k; $i++) {

	    my $ss = substr($s, $i, $k);
	    if ($ss !~ /\-/) {
		$H[$j]->{ $ss } += 1; 
	    }
	}

	
	#
	#  use a randomized version, ie shuffle both the initial and aligned sequence
	#

	if ($j == 1) {
	    my $sr = Sets::shuffle_seq($s);
	    my $l = length($sr);
	    for (my $i=0; $i<$l-$k; $i++) {
		my $ss = substr($sr, $i, $k);
		if ($ss !~ /\-/) {
		    $H[$j+1]->{ $ss } += 1; 
		}
	    }
	}

	
    }
}




foreach my $r (@$a_ref_kmers) {
    
    #print "HYP $r $H[2]->{$r}, $H[0]->{$r}, $H[1]->{$r}, N=\"$N\"\n\n";

    
    #my $todo = "hypergeom -i $H[2]->{$r} -s1 $H[0]->{$r} -s2 $H[1]->{$r} -N $N";
    #system($todo);
    my $f = $H[1]->{$r} / $H[0]->{$r};

    print "$r\t$H[0]->{$r}\t$H[1]->{$r}\t$f\t$H[2]->{$r}\n";
}
