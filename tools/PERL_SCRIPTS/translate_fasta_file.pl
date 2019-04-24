#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;
use Sets;
use strict;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my $i = 0;
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;

    my $st = Sets::translate($s);
    
    print ">$n\t$st\n\n";
    #foreach my $r (@$st) {
    #  print "$r => " . Sets::translate($r) . "\n";
    #}
}
