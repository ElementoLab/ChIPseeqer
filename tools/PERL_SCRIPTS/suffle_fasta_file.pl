# suffle all sequences in a FASTA file

use lib qw(/home/olly/PERL_MODULES);
use Sets;

use strict;

use Bio::SeqIO;




my $in  = Bio::SeqIO->new(-file => $ARGV[1] , '-format' => 'Fasta');

while ( my $o_seq = $in->next_seq() ) {
   
    my $seq_shuffled = Sets::shuffle($o_seq->seq);

    print ">";
    print $o_seq->display_id;
    print $seq_shuffled;
    print "\n\n";

}
  
