use lib qw(/home/olly/PERL_MODULES);

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    open OUT, ">tmp.seq";
    print OUT ">$n\n$s\n";
    close OUT;
    my $todo = "/home/olly/PERL_MODULES/PROGRAMS/blastz-source/blastz tmp.seq $ARGV[1] > blastz-$n.txt";
    system($todo);
}



