use lib qw(/home/olly/PERL_MODULES);

use Table;

my $ta = Table->new;

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my @OUT = sort { cmpk($b->[2], $a->[2]) } @$a_ref;

my $i = 1;
foreach my $r (@OUT) {
    $r->[0] = $i++;

    print join("\t", @$r) . "\n";
}


sub cmpk {
    my ($f1, $f2) = @_;

    if ($f1 =~ /inf/) {
	return 1;
    }

    if ($f2 =~ /inf/) {
	return -1;
    }
	
    return ($f1 <=> $f2);
    
}
