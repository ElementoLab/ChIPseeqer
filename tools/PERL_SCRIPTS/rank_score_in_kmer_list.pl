use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

my $a_ref = $ta->getArray();

my $i = 2;
foreach my $r (@$a_ref) {

    if ($r->[4] <= $ARGV[0]) {
	print "$i\n";
	exit();
    }
    
    $i ++;
}
