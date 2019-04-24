use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my @sums = ();
foreach my $r (@$a_ref) {
    
    my $n = scalar(@$r);
    
    for (my $i=0; $i<$n; $i++) {

	$sums[$i] += $r->[$i];
    }

}


print join("\t", @sums); print "\n";
