use lib qw(/home/elemento/PERL_MODULES);

use Fasta;
use Sets;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my @A = ();

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    print "$n\t$s\n";
    #push @A, $s;
}

