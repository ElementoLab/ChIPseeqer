use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;


my $ta = Table->new;

$ta->loadFile($ARGV[0]);
my $a_ref1 = $ta->getArray();



$ta->loadFile($ARGV[1]);
my $h_ref2 = $ta->getIndex(0);


foreach my $r (@$a_ref1) {

    if (defined($h_ref2->{$r->[0]})) {
	
	print "$r->[0]\t$r->[1]\t" . $h_ref2->{$r->[0]}->[1] . "\n";

    }

}
