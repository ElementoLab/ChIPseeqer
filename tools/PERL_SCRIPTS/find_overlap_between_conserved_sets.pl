use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;
use Hypergeom;


my $a_ref_f1 = Sets::getFiles($ARGV[0]);
my $a_ref_f2 = Sets::getFiles($ARGV[1]);

foreach my $f1 (@$a_ref_f1) {
    
    my $set1 = Sets::readSet($f1);

    foreach my $f2 (@$a_ref_f2) {

	my $set2 = Sets::readSet($f2);
	
	
	my $i  = Sets::getOverlapSize($set1, $set2);
	my $s1 = scalar(@$set1);
	my $s2 = scalar(@$set2);
	
	my $p = Hypergeom::cumhyper($i, $s1, $s2, 11265);


	if ($p < 1e-6) {
	    
	    
	    print Sets::basename($f1);
	
	    print "\t";
	    
	    print Sets::basename($f2);
	    
	    print "\t";
	    
	    print sprintf("%4.3e\n", $p);
	} 


    }

}
