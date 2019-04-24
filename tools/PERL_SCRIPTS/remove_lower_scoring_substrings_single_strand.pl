use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $n = scalar(@$a_ref);

my @a = ();
for (my $i=0; $i<$n-1; $i++) {
    for (my $j=$i+1; $j<$n; $j++) {
	my $ki = $a_ref->[$i]->[0];
	my $kj = $a_ref->[$j]->[0];

	#print "$kj =~ /$ki/ ?\n";
	
	if ($ki =~ /$kj/) {
	    $a[$j] = 1;
	} 
    
    }
}

for (my $i=0; $i<$n; $i++) {
    if (!defined($a[$i])) {
	print join("\t", @{$a_ref->[$i]}); print "\n";
    }

}
