#
#  input: 2x DNA sequence files
#
use lib qw(/home/elemento/PERL_MODULES);
use Fasta;
use ClustalW;

my $fa1 = Fasta->new;
$fa1->setFile($ARGV[0]);

my $fa2 = Fasta->new;
$fa2->setFile($ARGV[1]);

if ($ARGV[2]) {
    open OUT2, ">$ARGV[2]";
}

while (1)  {
    
    my $a_ref1 = $fa1->nextSeq();
    my $a_ref2 = $fa2->nextSeq();
    
    last if (!$a_ref1 || !$a_ref2);
    
    my ($n1, $s1) = @$a_ref1;
    my ($n2, $s2) = @$a_ref2;

    open OUT, ">fasta.seq";

    print OUT ">$n1.1\n$s1\n";
    print OUT ">$n2.2\n$s2\n";

    close OUT;

    my $todo = "clustalw fasta.seq -OUTORDER=INPUT -TYPE=DNA > /dev/null";
    system($todo);

    my $cl = ClustalW->new;
    $cl->setFile("fasta.aln");

    my $a_ref_aln = $cl->getInputOrderedSeqArray();

    my $n = length($a_ref_aln->[0]);

    my @a = split //, $a_ref_aln->[0];
    my @b = split //, $a_ref_aln->[1];

    print ">$n1\n";
    for (my $i=0; $i<$n; $i++) {
	if ($a[$i] eq $b[$i]) {
	    print $a[$i];
	} else {
	    print "-";
	}
    }
    print "\n\n";

    if ($ARGV[2]) {
	print OUT2 ">$n1\n$a_ref_aln->[0]\n\n";
    }
    

    #<STDIN>;
}

if ($ARGV[2]) {
    close OUT2; 
}

unlink "fasta.seq";
