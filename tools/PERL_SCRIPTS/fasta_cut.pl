#
# 1. load the lengths
#


use lib qw(/home/olly/PERL_MODULES);
use Table;
use Fasta;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0, 1);


my $deflen = $ARGV[2];


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
while ( my $a_ref = $fa->nextSeq ) {
    my ($name, $seq) = @{$a_ref};
    
    if (defined($h_ref->{$name})) {
	$seq = substr($seq, 0, $h_ref->{$name});
    } else {
	$seq = substr($seq, 0, $deflen);
    }

    print ">$name\n$seq\n\n";

}



