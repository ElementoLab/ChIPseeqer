use lib qw(/home/olly/PERL_MODULES);
use Table;


my $ta = Table->new;
$ta->setLimit(1000);
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
    
    my $nt1 = Sets::getNbLetterOccurrences($r->[0], 'T');
    my $nt2 = Sets::getNbLetterOccurrences($r->[1], 'T');

    if (($nt1 < length($r->[0])-1) && 
	($nt2 < length($r->[1])-1) ) {
	print join("\t", @$r); print "\n";
    }
}
