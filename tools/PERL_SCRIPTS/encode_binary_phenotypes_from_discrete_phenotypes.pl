use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r); print "\n";
my $i = 0;
my @names = ();
foreach my $r (@$a_ref) {
    
    my @a = ();
    my $H = ();
    foreach my $s (@$r) {
	
	if (defined($H{$s})) {
	   
	} else {
	    $H{$s} = $i++;
	    push @names, $s;
	    
	}

	push @a, $H{$s};
    }
    print join("\t", @a); print "\n";
}

print scalar(@names); print "\t";
print join("\t", @names); print "\n";


