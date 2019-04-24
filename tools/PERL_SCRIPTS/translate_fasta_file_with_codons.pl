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

    my $st = Sets::getArrayOfCodons($s);

    print ">$n\n";
    my $c = 1;
    foreach my $r (@$st) {
      print "$c: $r => " . Sets::translate($r) . "\n";
      $c++;
    }
    
}
