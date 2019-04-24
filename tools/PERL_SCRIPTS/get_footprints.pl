# input : gene, D mel region, get everything that aligned to that region
use lib qw(/home/olly/PERL_MODULES);
use Sets;

open IN, "OUT_FOOT/$ARGV[0].txt";
while (my $l = <IN>) {
    chomp $l;
    
    my @a = split /\t/, $l;

    if (Sets::sequencesOverlap($ARGV[1], $ARGV[2], $a[1], $a[2])) {
	
	# ok, there is an overlap between this footprint and the D mel region
	# now determine exactky what this overlap is in terms of sequences
	my $s = "";
	my @b = split //, $a[3];
	$cnt = $a[1];   # position in dmel
	for (my $i=0; $i<length($a[3]); $i++) {
	    if ($b[$i] ne '-') {

		if (($cnt >= $ARGV[1]) && ($cnt <= $ARGV[2])) {
		    $s .= $b[$i];
		}

		$cnt ++;

	       
	    }

	    
	}

	print ">F\n$s\n";
    }
    
}



