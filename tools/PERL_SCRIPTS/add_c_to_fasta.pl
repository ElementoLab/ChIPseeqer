use lib qw(/home/olly/PERL_MODULES);
use Fasta;


my $fa = Fasta->new;


$fa->setFile($ARGV[0]);

while ( my $a_ref = $fa->nextSeq ) {
    my ($name, $seq) = @{$a_ref};
    print ">C$name\n$seq\n\n";
}
