use lib qw(/home/olly/PERL_MODULES);

use Fasta;
use Sequence;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my %LENGTH  = ();
my %PROTEIN = ();

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;

    my $nc = $n;

    $nc =~ s/\-.+$//g;

    if (defined($LENGTH{$nc})) {
	if (length($s) > $LENGTH{$nc}) {
	    $LENGTH{$nc}  = length($s);
	    $PROTEIN{$nc} = $n;
	}
    } else {
	$LENGTH{$nc}  = length($s);
	$PROTEIN{$nc} = $n;

    }
}

my $se = Sequence->new;
$se->setBlastDB($ARGV[0]);
foreach my $n (keys(%PROTEIN)) {
    #print "$n\t$PROTEIN{$n}\n";

    my $seq = $se->getSequenceFromBlastDB($PROTEIN{$n}, 0, 0);

    print ">$n\n$seq\n\n";
}
