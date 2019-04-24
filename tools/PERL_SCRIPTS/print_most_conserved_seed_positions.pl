#
#
#  inout : mirnas, orthologous upstream regions
#  what : do a recompare of seed - gap - 3nt
#


use lib qw(/home/olly/PERL_MODULES);

use Fasta;
use Table;
use strict;

#
#  read in the kmers
#

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

my $h_ref = $ta->getIndex(0);
my $a_ref_kmers = $ta->getColumn(0);

#
#  read in the miRNAs
#
my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $cnt0 = 0;
my $cnt1 = 0;

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    my $l = length($s);

    my $ss = $s; $ss = uc($ss); $ss =~ s/U/T/g;

    print ">$n\t$ss\n";

    
    my @a_pos = ();
    

    foreach my $k (@$a_ref_kmers) {
	my $a_ref_pos = Sets::getREMotifPositions_singlestrand(Sets::getComplement($k), $ss);
	if (scalar(@$a_ref_pos) == 0) {
	    next;
	}
	my $pos    = shift @$a_ref_pos;  
	if ($pos > 1) {
	    next;
	}
	
       $a_pos[ $pos ] = $h_ref->{$k}->[4];
    }

    if ($a_pos[0] && $a_pos[1]) {
	if ($a_pos[0] > $a_pos[1]) {
	    $cnt0 ++;
	} else {
	    $cnt1 ++;
	}
    } elsif ($a_pos[0]) {
	$cnt0 ++;
    } elsif ($a_pos[1]) {
	$cnt1 ++;
    }
    
    
}

print "$cnt0\t$cnt1\n";
