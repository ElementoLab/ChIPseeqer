#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
 
    my @ss = split /N+/, $s;

    my $i = 1;
    foreach my $m (@ss) {
      print ">$n-$i\n$m\n";
      $i++;
    }
    
}
