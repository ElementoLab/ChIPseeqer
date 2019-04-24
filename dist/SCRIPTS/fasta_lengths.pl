#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while ( my $a_ref = $fa->nextSeq ) {
    my ($name, $seq) = @{$a_ref};
 
    my $l = length($seq);
    
    print "$name\t$l\n";
}
