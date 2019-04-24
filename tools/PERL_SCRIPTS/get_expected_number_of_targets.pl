use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;

my $a_ref_files = Sets::readSet($ARGV[0]);

my @a1 = ();
my @a2 = ();
foreach my $f (@$a_ref_files) {
    
    $ta->loadFile($f);
    my $a_ref = $ta->getArray();

    my $i =0;
    foreach my $r (@$a_ref) {
	
	#print "$r->[0]\n" if ($i == 0);
	
	$a1[$i] += $r->[0];
	$a2[$i] = $r->[1];
	$i++;
    }
    


}

my $nf = scalar(@$a_ref_files);
my $nr = scalar(@a1);


for (my $i=0; $i<$nr; $i++) {
    
    $a1[$i] = $a1[$i] / $nf;

    print sprintf("%3.2f", $a1[$i]) . "\t$a2[$i]\n";
    
}
