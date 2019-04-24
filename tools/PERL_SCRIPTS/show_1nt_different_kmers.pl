use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();


my $k = $ARGV[0];


for (my $i=0; $i<length($k); $i++) {
    
    my $pat = $k; my @a = split //, $pat; 
    $a[$i]='.'; 
    $pat    = join('', @a);
    my $j = 1;
    foreach my $r (@$a_ref) {

        if ($r->[0] eq $k) {
	    $j++; next;
	}

	if ($r->[0] =~ /$pat/) {
	    print "$k -> $r->[0]\t$j\n";

	} 
	
	$j++;
    }

}
