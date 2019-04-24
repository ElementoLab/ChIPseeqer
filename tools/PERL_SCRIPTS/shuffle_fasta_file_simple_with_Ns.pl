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

    my $ss = $s;

    $s =~ s/N//g;

    my $seq_shuffled = Sets::shuffle_seq($s);

    my @a1 = split //, $ss;
    my @a2 = split //, $seq_shuffled;

    my $news = "";
    my $j = 0;
    for (my $i=0; $i<@a1; $i++) {
      if ($a1[$i] eq 'N') {
	$news .= 'N';
      } else {
	$news .= $a2[$j];
	$j ++;
      }
    }

    #print ">";
    #print $o_seq->display_id;
    #print "\n";
    
    print ">$n\n$news\n\n";
    #print "\n\n";

}
