use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;

my $a_ref_files = Sets::readSet($ARGV[0]);

die "please provide a directory\n" if (!$ARGV[1]);

my @a1 = ();
my @a2 = ();
foreach my $f (@$a_ref_files) {
    
    $ta->loadFile("$ARGV[1]/$f");
    my $a_ref = $ta->getArray();

    my $i =0;
    foreach my $r (@$a_ref) {
	
	#print "$r->[0]\n" if ($i == 0);
	
	push @{ $a1[$i] }, $r->[0];
	$a2[$i] = $r->[1];
	$i++;
    }
    


}

my $nf = scalar(@$a_ref_files);
my $nr = scalar(@a1);


for (my $i=0; $i<$nr; $i++) {
    
    #$a1[$i] = $a1[$i] / $nf;

    print "$a2[$i]\t";

    my $m = Sets::average($a1[$i]);
    my $s = Sets::stddev($a1[$i]);
    
    print "$m\t$s\n";

    #print join("\t", @{ $a1[$i] }); print "\n";


    #print sprintf("%3.2f", $a1[$i]) . "\t$a2[$i]\n";
    
}
