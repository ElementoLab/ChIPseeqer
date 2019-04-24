#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);
use Sets;
use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $r = $fa->nextSeq()) {
    my ($n, $s) = @$r;
    my $sc = undef;
    if (($s =~ /U/) && ($s !~ /T/)) {
      $sc = Sets::getRNAComplement($s);
    } else {
      $sc = Sets::getComplement($s);
    }
    
    print ">$n\n$sc\n\n";

}






 
