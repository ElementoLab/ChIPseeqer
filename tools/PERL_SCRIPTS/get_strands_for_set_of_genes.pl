use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

my $h_ref = $ta->getIndex(4);

my $a_ref = Sets::readSet($ARGV[0]);

my $a = 0;
my $b = 0;
foreach my $r (@$a_ref) {
    
    print "$r\t" . $h_ref->{$r}->[3] . "\n";
    
    if ($h_ref->{$r}->[3] == 1) {
	$a++;
    } else {
	$b++;
    }
}


print "$a\t$b\n";
