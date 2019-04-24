use lib qw(/home/olly/PERL_MODULES);
use Bio::SeqIO;
use Bio::Seq;
use Sequence;



my $seqio  = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'genbank');
my @feats = $seq->all_SeqFeatures;


my @CLUSTER = ();


my $s = Sequence->new;
$s->setBlastPath("/home/olly/COMPARATIVE_YEAST/ORFS/BLAST");
$s->setBlastDB($ARGV[1]);

# temporary storage files
my $s_tmpreport = "/tmp/fasta.report";
my $s_tmpstore  = "/tmp/fasta.store";

my $s_fastapath = "/home/olly/MPLPC18/STUDENTS/ISABELLE/FASTA";

foreach my $feat (@feats) {
    my $tag = $feat->primary_tag();
    if ($tag eq "gene") {

	# cree un nouveau gène
	my %GENE = ("GENE"  => join ('', $feat->each_tag_value("gene")),
		    "START" => $feat->start,
		    "END"   => $feat->end);
	
	# get the sequence for this gene 
	my $seq = $s->getSequenceFromBlastDB($ARGV[2], $GENE{START}, $GENE{END});
		
	#
	# compare this sequence to all the other ones in the same cluster
	#

	# store this sequence
	open SEQ, ">$s_tmpstore1";
	print SEQ ">SEQ1\n$seq\n\n";
	close SEQ;
	
	# flag to decide if we store or not
	my $flag = 0;

	foreach my $c (@CLUSTER) {

	    # store this sequence
	    open SEQ, ">$s_tmpstore2";
	    print SEQ ">SEQ2\n" . $c->{SEQ} . "\n\n";
	    close SEQ;
	    
	    # use FASTA to align the current sequence against the reference sequence 
	    my $output = `$s_fastapath/fasta33 -A -Q -b1 $s_tmpstore1 $s_tmpstore2`;
	    
	    # output the report
	    open OUT, ">$s_tmpreport";
	    print OUT $output;
	    close OUT;
	    
	    
	    print $s_tmpreport;
	    
	    # create a BioPerl object to handle the report
	    my $searchio = new Bio::SearchIO(-format => 'fasta',
					     -file   => $s_tmpreport);
	    my $result = $searchio->next_result();
	    my $hit    = $result->next_hit;
	    my $hsp = $hit->next_hsp;
	    
	    my $ident  = $hsp->frac_identical;
	    
	    # get the query string
	    my $s_seq = $hsp->query_string;
	    
	    # query string length / initial query string length
	    my $fraclen = length($s_seq) / length($c->{SEQ});
	    
	    
	    if ( ($ident > 0.7) && ($fraclen > 0.6) ) {
		# add to the cluster
		$flag = 1;
	    } else {
		$flag = 0;
		
	    }
	    
	}
	
	if ($flag == 1) {
	    # first store the sequence 
	    $GENE->{SEQ} = $seq;

	    # then add to the cluster
	    push @CLUSTER, \@GENE;
	} else {
	    # stop, new cluster
	    print "Current cluster has " . scalar(@CLUSTER). " tandem genes\n";
	    @CLUSTER = ();
	}
    }

    
}

