use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;

# load kmers
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $n = scalar(@$a_ref);

my $w = (defined($ARGV[1])?$ARGV[1]:100);

for (my $i=0; $i<$n-$w; $i++) {
    
    my $cnt = 0;
    for (my $j=$i; $j<$i+$w; $j++) {
	
	$cnt += $a_ref->[$j]->[1];
	
    }

    $cnt = $cnt / $w;


    print $a_ref->[$i]->[0]; print "\t";
    print "$cnt\n";
    
}
