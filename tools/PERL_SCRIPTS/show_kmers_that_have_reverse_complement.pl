use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $h_ref = $ta->getIndex(0);
my @keys  = keys(%$h_ref);

foreach my $r (@keys) {
    
    
    my $c = Sets::getComplement($r);

    if (defined($h_ref->{ $c })) {
	
	print "$r\n";
    }
    
}

