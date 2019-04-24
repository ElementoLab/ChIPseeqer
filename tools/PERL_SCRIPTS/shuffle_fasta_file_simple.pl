# suffle all sequences in a FASTA file
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;
use strict;

use Fasta;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;

    my $seq_shuffled = Sets::shuffle_seq($s);

    #print ">";
    #print $o_seq->display_id;
    #print "\n";
    
    print ">$n-shu\n$seq_shuffled\n\n";
    #print "\n\n";

}
