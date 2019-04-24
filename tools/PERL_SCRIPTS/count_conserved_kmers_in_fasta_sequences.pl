use lib qw(/home/olly/PERL_MODULES);

use Fasta;
use Sets;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my @A = ();

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    push @A, $s;
}


my $s1 = shift @A;
my $l  = length($s1);
my $n  = scalar(@A);

for (my $i=0; $i<$l-6; $i++) {
    
    my $kmer = substr($s1, $i, 6);

    next if ($kmer =~ /\-/);
    
    my $cnt = 0;
    foreach my $s (@A) {
	
	$okmer = substr($s, $i, 6);
	
	$cnt++ if ($okmer eq $kmer);

    }

    if ($cnt >= 3) {
	
	my $mi = Sets::isMiRNATarget($kmer, "fly");
	
	my $txt = join('/', @$mi);

	print "$kmer\t$cnt\t$txt\n";
    } 
    
}

