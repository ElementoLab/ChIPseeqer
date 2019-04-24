use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;


my $ta1 = Table->new;

$ta1->loadFile($ARGV[0]);
my $n1 = $ta1->getNbColumns();

my $ta2 = Table->new;

$ta2->loadFile($ARGV[1]);
my $n2 = $ta2->getNbColumns();

for (my $i=1; $i<$n1; $i++) {
    
    
    my $c1_ref = $ta1->getColumn($i);
    my $l1     = shift @$c1_ref;

    for (my $j=1; $j<$n2; $j++) {
	
	my $c2_ref = $ta2->getColumn($j);
	my $l2     = shift @$c2_ref;
	
	my $c      = Sets::pearson($c1_ref, $c2_ref);
    
	if (abs($c) >= 0.1) {
	    print "$l1\t$l2\t$c\n";
	}
    }
    
    
    
    
}
