use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;
use strict;

my $t = Table->new;
$t->setLimit(800);
$t->loadFile($ARGV[0]);

my $a_ref = $t->getArray;


my $n  = scalar(@$a_ref);

for (my $i=0; $i<$n; $i++) {
    #print "$i: $a_ref->[$i]->[0]\n";
    for (my $j=$i+1; $j<$n; $j++) {
	
	if ($a_ref->[$i]->[0] eq Sets::getComplement($a_ref->[$j]->[0])) {
	    $a_ref->[$j]->[5] = 'N';
	
	    #print "$j:Eliminate $a_ref->[$j]->[0] ..\n";
	}

    }

}

foreach my $r (@$a_ref) {
    if ($r->[5] ne 'N') {
	print join("\t", @$r); print "\n";
    }
}
