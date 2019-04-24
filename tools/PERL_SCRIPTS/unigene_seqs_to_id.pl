use lib qw(/home/olly/PERL_MODULES);
use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;

    my ($hs) = $n =~ /ug=(.+?)\ /;

    print ">$hs\n$s\n\n";
}

