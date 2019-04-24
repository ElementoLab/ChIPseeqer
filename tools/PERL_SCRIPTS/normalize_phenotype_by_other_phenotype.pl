#
#  different phenotypes
#

use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;

$ta->loadFile($ARGV[0]);
my $a_ref1 = $ta->getArray();

$ta->loadFile($ARGV[1]);
my $h_ref2 = $ta->getIndex(0);

my $a_ref2   = $h_ref2->{"Nb_genes"};
shift @$a_ref2;

my $firstline = shift @$a_ref1;

print join("\t", @$firstline); print "\n";

foreach my $r (@$a_ref1) {
    
    print shift(@$r); print "\t";
    my @row = ();
    for (my $i=0; $i<scalar(@$r); $i++) {
	if (($r->[$i] eq "inf") || ($a_ref2->[$i] eq "inf")) {
	    push @row, "inf";
	} else {
	    push @row, sprintf("%4.3f", $r->[$i]/$a_ref2->[$i]);
	}
    }
    print join("\t", @row); print "\n";

}



