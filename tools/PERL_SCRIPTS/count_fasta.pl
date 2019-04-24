#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);
use Fasta;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
my $i = 0;
while ( my $a_ref = $fa->nextSeq ) {
    $i++;
}


print "$i\n";
