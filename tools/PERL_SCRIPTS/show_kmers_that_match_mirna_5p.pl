use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;

# load kmers
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


# load 5' kmers
my $h_ref_5p = Sets::getIndex($ARGV[1]);

foreach my $r (@$a_ref) {
    
    if ( defined($h_ref_5p->{$r->[0]}) ) {
	print "$r->[0]\t1\n";
    } else {	
	print "$r->[0]\t0\n";
    }
    
}
