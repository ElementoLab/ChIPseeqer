use lib qw(/home/olly/PERL_MODULES);
use Markov;

my $f   = $ARGV[0];
my $len = $ARGV[1];

open IN, $f;

while (my $l = <IN>) {


    chomp $l;

    print "$l\n" if ($l =~ /^\>/);
    
    $seq .= $l;

    if (length($seq) >= $len) { 
	my $o_markov = Markov->new();
	$o_markov->calcFrequenciesFromSeq($seq);
	my $newseq = $o_markov->generate1rstOrder(length($seq));
	print "$newseq\n";
	
	$seq = "";
    } 


}
my $o_markov = Markov->new();
$o_markov->calcFrequenciesFromSeq($seq);
my $newseq = $o_markov->generate1rstOrder(length($seq));
print "$newseq\n";
close IN;
