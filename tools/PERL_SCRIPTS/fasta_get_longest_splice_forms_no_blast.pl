BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Fasta;
use Sequence;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my %LENGTH  = ();
my %PROTEIN = ();
my %SEQUENCE = ();

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;

    my $nc = $n;

    $nc =~ s/\-.+$//g;

    if (defined($LENGTH{$nc})) {
	if (length($s) > $LENGTH{$nc}) {
	    $LENGTH{$nc}  = length($s);
	    $PROTEIN{$nc} = $n;
	    $SEQUENCE{$nc} = $s;
	}
    } else {
	$LENGTH{$nc}  = length($s);
	$PROTEIN{$nc} = $n;
	$SEQUENCE{$nc} = $s;

    }
}


foreach my $n (keys(%PROTEIN)) {
    
    print ">$n\n$SEQUENCE{$n}\n\n";
}
