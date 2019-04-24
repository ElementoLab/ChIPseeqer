use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;
use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $nb     = $ARGV[1];
my $si     = $ARGV[2];
my $suffix = $ARGV[3];

my @SEQS = ();
while (my $a_ref = $fa->nextSeq()) {

    my ($n, $s) = @$a_ref;
    
    push @SEQS, $s;
    
}



srand;

for (my $i=1; $i<=$nb; $i++) {

    my $r = int(rand($nb));

    my $l = int(rand(length($SEQS[$r])-$si));

    my $s = substr($SEQS[$r], $l, $si) . $suffix;

    print ">rand$i\n$s\n";

}
