use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


#print "sorting ...\n";

my @a_sorted = sort { $a->[3] <=> $b->[3] } @$a_ref;

#print "sorted.\n";

my $n = scalar(@a_sorted);


my @DEAD = ();

for (my $i=0; $i<$n; $i++) {

    next if ($DEAD[$i] == 1);

    
    
    my $t = $a_sorted[$i]->[2];
    my $s = $a_sorted[$i]->[3];
    my $e = $a_sorted[$i]->[4];
    my $r = $a_sorted[$i]->[6];
    #print "looking at PRE-TRANSPOSON $s -> $e\n";

    for (my $j=$i+1; $j<$n; $j++) {

	#print "  does  $a_sorted[$j]->[3], $a_sorted[$j]->[4] overlap with $s, $e? ";
	
	if (Sets::sequencesOverlap($s, $e, $a_sorted[$j]->[3], $a_sorted[$j]->[4])) {

	    if (($t eq $a_sorted[$j]->[2]) && ($r eq $a_sorted[$j]->[6])) {

		$s = Sets::min($s, $a_sorted[$j]->[3]);
		$e = Sets::max($e, $a_sorted[$j]->[4]);
		#print "YES\n";
		$DEAD[$j] = 1;
		
	    } else {
		#print "YES, but different type $a_sorted[$j]->[2] or strand $a_sorted[$j]->[6]\n";
	    }

	} else {
	    #print "NO\n";
	    print "$t\t$s\t$e\n";

	    last;
	}
	
	
    }

}
