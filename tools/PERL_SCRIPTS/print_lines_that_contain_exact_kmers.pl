use lib qw(/home/olly/PERL_MODULES);

use Table;
use strict;

# load kmers
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref_kmers = $ta->getIndex(0);

# load file
$ta->loadFile($ARGV[1]);
my $a_ref_file = $ta->getArray();

my $cnt = 1;
my $renum = $ARGV[2];

foreach my $r (@$a_ref_file) {
    foreach my $s (@$r) {
	my $ss = $s;
	$ss =~ s/\ //g;

	

	if ($h_ref_kmers->{$ss}) {
	    
	    if ($renum == 1) {
	    $r->[0] = $cnt++;
	}

	    
	    print join("\t", @$r); print "\n";
	}
    }
}
