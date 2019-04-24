BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while ( my $a_ref = $fa->nextSeq ) {
    my ($name, $seq) = @{$a_ref};
 
    my $subseq = substr($seq, length($seq)-$ARGV[1], $ARGV[1]);
    
    print ">$name\n$subseq\n\n";
}
