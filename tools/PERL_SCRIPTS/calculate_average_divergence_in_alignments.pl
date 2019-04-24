use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my @D = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    my $l = length($s);

    my @a = split //, $s, -1;
    my $cnt = 0;
    foreach my $i (@a) {
	if ($i ne "-") {
	    $cnt ++;
	}
    }

    push @D, $cnt/$l;
}



print Sets::average(\@D); print "\n";
