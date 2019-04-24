use lib qw(/home/olly/PERL_MODULES);
use Table;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {
    if (!defined($H{$r->[0]})) { 
	print join(",", @$r); print "\n";
	$H{$r->[0]} = 1;
    } else {
	#$H{$r->[0]} = 1;
    }
}
