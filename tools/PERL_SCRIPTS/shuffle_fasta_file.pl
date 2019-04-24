# suffle all sequences in a FASTA file

use lib qw(/home/olly/PERL_MODULES);
use Bio::SeqIO;
use Sets;
use strict;
use Markov;

use Fasta;

my $in  = Bio::SeqIO->new(-file => $ARGV[0] , '-format' => 'Fasta');

while ( my $o_seq = $in->next_seq() ) {
   
    my $o_markov = Markov->new();
    
    if ($o_seq->seq !~ /[ACTG]/) {
	print ">" . $o_seq->id() . "\n" . $o_seq->seq . "\n\n"; 
    } else {
	$o_markov->calcFrequenciesFromSeq($o_seq->seq);
	my $s_new = $o_markov->generate1rstOrder($o_seq->length);
	print ">" . $o_seq->id() . "\n" . $s_new . "\n\n"; 
    } 

    #my $seq_shuffled = Sets::shuffle_seq($o_seq->seq);

    #print ">";
    #print $o_seq->display_id;
    #print "\n";
    
    #print $seq_shuffled;
    #print "\n\n";

}
  
