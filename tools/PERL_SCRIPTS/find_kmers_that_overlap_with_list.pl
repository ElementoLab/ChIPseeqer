use lib qw(/home/olly/PERL_MODULES);
use Table;

my $a_ref_k = Sets::readSet($ARGV[0]);


my %IDX = ();
foreach my $k (@$a_ref_k) {
    $IDX{ $k } = 1;

}
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();


foreach my $k (@$a_ref_k) {
    
    

    foreach my $r (@$a_ref) {
	
	next if ($k eq $r->[0]);
	next if ($IDX{ $r->[0] });
	
	if (Sets::seqOverlap($k, $r->[0]) == 6) {
	    print "$k\t$r->[0]\n";
	}
	

    }
    
}
