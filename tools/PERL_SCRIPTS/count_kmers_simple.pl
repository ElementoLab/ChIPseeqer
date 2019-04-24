#
#  input: miRNAs
#
#

use lib qw(/home/olly/PERL_MODULES);
use Sets;

use Fasta;
use strict;
use Hypergeom;

my $a_ref_kmers = Sets::readSet($ARGV[0]);
my $const = 500;
my @H = ();
my %CO = ();

my $k = undef;
foreach my $r (@$a_ref_kmers) {
    $H[0]->{ $r } = 0;
    $H[1]->{ $r } = 0;
    $H[2]->{ $r } = 0;
    $CO{ $r } = 0;
    $k = length($r);
}


my $N = 0;
my %LEN = ();
for (my $j=0; $j<3; $j++) {
    my $fa = Fasta->new;
    $fa->setFile($ARGV[$j+1]);
    
    $N  = 0;

    while (my $a_ref = $fa->nextSeq()) {
	my ($n, $s) = @$a_ref; $N++;
	my $l = length($s); if ($j == 0) { $LEN{ $n } = $l };
	my %L = ();
	for (my $i=0; $i<$l-$k; $i++) {
	    
	    my $ss = substr($s, $i, $k);
	    
	    
	    if (($ss !~ /\-/) && (!defined($L{$ss}))) {
		$H[$j]->{ $ss } += 1; 

		if ($j == 2) {
		    $CO{ $ss } ++;
		}

		$L{ $ss } = 1;
	    }
	    
	}
	
    }
}




foreach my $r (@$a_ref_kmers) {
    
    #print "HYP $r $H[2]->{$r}, $H[0]->{$r}, $H[1]->{$r}, N=\"$N\"\n\n";

    
    #my $todo = "hypergeom -i $H[2]->{$r} -s1 $H[0]->{$r} -s2 $H[1]->{$r} -N $N";
    #system($todo);
    my $p = -1 * Hypergeom::lcumhyper($H[2]->{$r}, $H[0]->{$r}, $H[1]->{$r}, $N);


    print "$r\t$H[0]->{$r}\t$H[1]->{$r}\t$H[2]->{$r}\t$p\t$CO{$r}\n";
}
