use lib qw(/home/olly/PERL_MODULES);
use Table;

my $f1 = shift @ARGV;
my $ta = Table->new;

my %H = ();
foreach my $f (@ARGV) {
    $ta->loadFile($f);
    my $a_ref = $ta->getArray();

    foreach my $r (@$a_ref) {
	$H{ $r->[0] } = $r;
    }
    
    
}

$ta->loadFile($f1);
my $a_ref = $ta->getArray();


foreach my $r (@$a_ref) {
        
    my $c = Sets::getComplement($r->[0]);

    print "$r->[0]\t$r->[1]\t$r->[2]\t" . int(0.5+$r->[4]) . "\t$H{$c}->[0]\t$H{$c}->[1]\t$H{$c}->[2]\t" . int(0.5+$H{$c}->[4]) . "\n";
    
}

