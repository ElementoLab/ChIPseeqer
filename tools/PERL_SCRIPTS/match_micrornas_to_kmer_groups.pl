use lib qw(/home/olly/PERL_MODULES);

use Sets;
use Table;
use Fasta;
use strict;

if (!$ARGV[1]) {
	die "Usage : prg seq kmers\n";

}

my @mirnas = ();
my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
while (my $a_seq = $fa->nextSeq()) {
    my @a = @$a_seq;
    push @mirnas, \@a;
}

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my $only5       = $ARGV[2];
my $fulldisplay = 1;

my $i = 1;
foreach my $r (@$a_ref) {
    
    #my @g = split /\t/, $l, -1;

    print join("\t", @$r); print "\n"; 
    
    my $match = 0;
    my $rank  = 1;

    foreach my $m (@mirnas) {
	my ($n, $s) = @$m;
	$s =~ s/t/u/g;
	my $ss = uc($s);
	$ss =~ s/U/T/g;
	$ss = Sets::getComplement($ss);
	
	#
	# traverse all k-mers
	#
	foreach my $k0 (@$r) {
	    
	    my $k = $k0;

	    my $kmer =  substr($k, 0, length($k));
	    my $km   =  $kmer;
	    $km      =~ s/N/\./g;
	    
	    #
	    # does the mirna contain the kmer 
	    #
	    if ($ss =~ /$km/) {
		
		$km     =  Sets::getComplement($km);
		$km     =~ s/N/\./g;
		my $lkm =  lc($km);
		$lkm    =~ s/t/u/g;
		my $sc = $s;
		$sc      =~ s/$lkm/$km/g;
		my $a_ref_pos = Sets::getREMotifPositions($km, $sc);
		my $p   = shift @$a_ref_pos;
		if (defined($only5) && ($p > 2)) {
		    $rank++;
		    next;
		}
		
		print "\t$k\t$rank:$n\t$sc\t$p\n"; 

		$match++ if ($p <= 2);
	    }
	    
	   
	    $rank ++;
	    
	    #last if ($rank > 600);
	}
	
	#$fa->dispose();
	
	if ($fulldisplay == 0) {
	    $r->[5] = $match;
	    print join("\t", @$r); print "\n";
	}
	
	
	$i++;
    }
}
    
