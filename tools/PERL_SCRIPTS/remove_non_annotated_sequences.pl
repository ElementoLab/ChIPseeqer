#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);
use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);


while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    
    next if (($seq =~ /unavailable/) || ($seq =~ /annotated/));
    
    
    print ">$name\n$seq\n\n";
    
}
