use lib qw(/home/olly/PERL_MODULES);

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref1 = $ta->getIndex(0);


$ta->loadFile($ARGV[1]);
my $a_ref2 = $ta->getArray;


#my $cnt = 1;
foreach my $r (@$a_ref2) {
    #$r->[1] =~ s/\ //g;

    if ($h_ref1->{ $r->[0] }) {

	#$r->[0] = $cnt; $r->[1] = sprintf("%10s", $r->[1]);
	print join("\t", @$r) . "\n";

	#$cnt++;
    }

    
}
