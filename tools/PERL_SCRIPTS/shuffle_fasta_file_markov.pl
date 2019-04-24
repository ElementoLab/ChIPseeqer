#!/usr/bin/perl

# suffle all sequences in a FASTA file

use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;
use Markov;
use Fasta;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;

    
    #if (length($s) < 
       
    my $o_markov = Markov->new();
    
    #$o_markov->useN(1);
    
    $o_markov->calcFrequenciesFromSeq($s);
    my $s_new = $o_markov->generate1rstOrder(length($s));

    print ">$n-shu\n$s_new\n\n";

    #my $seq_shuffled = Sets::shuffle_seq($s);

    #print ">";
    #print $o_seq->display_id;
    #print "\n";
    
    #print ">$n\n$seq_shuffled\n\n";
    #print "\n\n";

}
