use lib qw(/home/olly/PERL_MODULES);
use MyBlast;
use Fasta;
use Sets;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

system(">$ARGV[0].rev");

my $tmpfile1 = "/tmp/fasta.1";
my $tmpfile2 = "/tmp/fasta.2";

while (my $a_ref = $fa->nextSeq) {
    
    my ($name, $seq) = @$a_ref;
    
    $fa->writeSeq($tmpfile1, $name, $seq);

    system("cat $tmpfile1 $ARGV[0].rev > $tmpfile2");

    system("mv $tmpfile2 $ARGV[0].rev");

    
}
