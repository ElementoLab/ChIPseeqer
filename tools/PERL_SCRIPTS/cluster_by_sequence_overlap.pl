use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;

my $ta = Table->new;
if (defined($ARGV[1])) {
  $ta->setLimit($ARGV[1]);
}
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $n = scalar(@$a_ref);
my $r = undef;

for (my $i=0; $i<$n; $i++) {

    next if ($taken[ $i] == 1);

    print "# $a_ref->[$i]->[0]\n";
    my %H = ();

    while (1) {

	$change = 0;
	for (my $j=0; $j<$n; $j++) {
	    
	    next if ($i == $j);
	    next if ($taken[ $j] == 1);
	 
	    #print "  compare $a_ref->[$i]->[0], $a_ref->[$j]->[0]\n";
   
	    
	    if (Sets::seqOverlap($a_ref->[$i]->[0], $a_ref->[$j]->[0], \$r) >= Sets::min(length($a_ref->[$i]->[0]), length($a_ref->[$j]->[0]))-1) {
		#print "  overlap with $a_ref->[$j]->[0]\n";
		$a_ref->[$i]->[0] = Sets::mergeOverlappingSequences($a_ref->[$i]->[0], $a_ref->[$j]->[0], $r);
		print "  --> extend $a_ref->[$i]->[0]\n";
		$taken[ $j ] = 1;
		$change = 1;
		last;
	    } else {
		
		my $c = Sets::getComplement($a_ref->[$j]->[0]);
		

		if (Sets::seqOverlap($a_ref->[$i]->[0], $c, \$r) >= Sets::min(length($a_ref->[$i]->[0]), length($a_ref->[$j]->[0]))-1) {
		    #print "  --> $c (overlap) $r\n";
		    $a_ref->[$i]->[0] = Sets::mergeOverlappingSequences($a_ref->[$i]->[0], $c, $r);
		    print "  --> extend $a_ref->[$i]->[0]\n";
		    $taken[ $j ] = 1;
		    $change = 1;
		    last;
		}
	    }
	    my $ov = Sets::getBestOverlapDP($a_ref->[$i]->[0], $a_ref->[$j]->[0], \$strand, \$decal, 1, 0, 1);

	    #print "     ov $a_ref->[$i]->[0], $a_ref->[$j]->[0] =$ov\n";
	    
	    if ($ov >= Sets::min(length($a_ref->[$i]->[0]), length($a_ref->[$j]->[0]))-1) {
		print "  --> add $a_ref->[$j]->[0]\n";
		$taken[ $j ] = 1;
		#$change = 0;
	    }
	}

	last if ($change == 0);
    }
}
