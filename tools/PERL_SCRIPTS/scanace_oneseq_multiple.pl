use lib qw(/home/olly/PERL_MODULES);
require "ScanACE.pm";
require 'Sequence.pm';
require 'Yeast.pm';

my $y = Yeast->new;

$list = "/home/olly/listing_simple.txt";

$db[0] = "/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/DATABASES/UTR5/utr5_sc_1000.fasta";
#$db[3] =  "/home/olly/YEAST_GENOMES/utr5_1000_mikatae_mit.fasta";
#$db[2] =  "/home/olly/YEAST_GENOMES/utr5_1000_bayanus_mit.fasta";
#$db[1] =  "/home/olly/YEAST_GENOMES/utr5_1000_paradoxus_mit.fasta";

#$db[2] =  "/home/olly/YEAST_GENOMES/utr5_1000_bayanus_mit.fasta";
#$db[2] =  "/home/olly/YEAST_GENOMES/utr5_1000_bayanus_mit.fasta";

open IN, $list or die "cannot find $list\n";
@a_genes = <IN>;
chomp @a_genes;
close IN;



my $o_seq = Sequence->new;
$o_seq->setBlastPath("/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/PROGRAMS/BLAST");

foreach my $s (@a_genes) {

    #print "$s\n";
    
    my @a = split /\t/, $s;

    next if (!$a[0]);

    print $a[0];
    print "\t";
    
    $ref = $y->getORF($a[0]);
    
    print $ref->{GENENAME}; #$y->getGeneName($a[0]);
    print "\t";
    print $ref->{ESSENTIAL};

    $chip = $y->getCHIPdata($ARGV[1], $a[0]);
    print "\t";
    print $chip->{RATIO};

    my $i = 0;
    foreach $g (@a) {
	
	
    
	$o_seq->setBlastDB($db[$i]);
	
	$i ++;

	next if (!$g);

	#$i ++;

	my $s1 = $o_seq->getSequenceFromBlastDB($a[0], 0, 0);
	
	if ($s1 && (length($s1) > 100)) {
	    open  OUT, ">tmp";
	    print OUT  ">$a[0]\n$s1\n\n";
	    close OUT;
	    
	} else {
	    
		#$i++;
	    next;
		#last;
	} 


	
	$s = ScanACE->new;
	
#$s->readMotif("../MOTIFS/RAP1.motif");
	$s->setNbMotifs(100);
#$s->setStdDev(3.0);
	$s->setGC(0.33);
#$s->setRun(0);
#$s->setThresholdScore(10.0);
	
	$s->runScanACE("tmp", $ARGV[0]);
	

#my $cnt = $s->getNbMotifsInFile($ARGV[0]);

	my $a_ref = $s->getSites;

	$top = 0.0;
	
	$tot = 0.0;    
	foreach my $a (@$a_ref) {
	    
	    #print "$ARGV[0]\t$a->[0]\t$a->[3]\t$a->[4]\t$a->[1]\t$a->[2]\n";
	    
	    $tot += $a->[4];

	    if ($a->[4] > $top) {
		$top = $a->[4];
	    } 
	}
	
	
	print "\t$tot\t$top";
	last;
    }

    print "\n";
}
