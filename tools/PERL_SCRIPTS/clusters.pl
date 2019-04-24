use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use Bio::Tools::pSW;
use Bio::Seq;
use Bio::Factory::EMBOSS;

GetOptions ('chr=s' => \$chr,
	    'dist=s'   => \$dist);



$seqio  = Bio::SeqIO->new(-file => $chr, '-format' => 'genbank');

$seq = $seqio->next_seq();

@feats = $seq->all_SeqFeatures;

# liste des clusters
@clusters = ();

# liste des gènes dans le cluster
my @clust = ();


# cluster courant 

$last_end = 0;

foreach $feat (@feats) {
    $tag = $feat->primary_tag();
    if ($tag eq "gene") {
	#print join ('', $feat->each_tag_value("gene")) . " (" . $feat->start . ", " .  $feat->end . ")\n";
	# cree un nouveau gène
	my %curclust = ("GENE"  => join ('', $feat->each_tag_value("gene")),
			"START" => $feat->start,
			"END"   => $feat->end,
			"SEQ"   => $feat->seq->seq);
	
	# si la distance entre le debut de ce gène et la fin du dernier gène <= delta, 
	if (($feat->start - $last_end) <= $dist) {
	    push @clust, \%curclust;
	    # sinon cree un nouveau cluster
	} else {
	    
	    my @copy =  @clust;
	    push @clusters, \@copy;
	    
	    @clust = ();
	    push @clust, \%curclust;
	}
	    
	$last_end = $feat->end;
    } 

    
}


# liste les clusters
foreach $c (@clusters) {
    @clust = @$c;
    print "Start new cluster\n";
    foreach $g (@clust) {
	print $g->{GENE} . " ";
    }

    print "\n";
}
    

$factory = new Bio::Tools::pSW( '-matrix' => 'blosum62.bla',
				'-gap' => 15,
				'-ext' => 5, );


print "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";

# liste les clusters
foreach $c (@clusters) {
    @clust = @$c;
    print "Start new cluster\n";
    for ($i=0; $i<scalar(@clust); $i++) {
	$o_seq1   = Bio::Seq->new( -seq => $clust[$i]->{SEQ});
	$o_aaseq1 = $o_seq1->translate();
	for ($j=$i+1; $j<scalar(@clust); $j++) {
	    $o_seq2   = Bio::Seq->new( -seq => $clust[$j]->{SEQ});
	    $o_aaseq2 = $o_seq2->translate();
	    
	    #print "Compare " . $clust[$i]->{SEQ} . " ***and*** " . $clust[$j]->{SEQ} . "..\n";
	    #print "Compare " . $o_seq1->seq() . " ***and*** " . $o_seq2->seq() . "...\n";
	    #$factory->align_and_show($o_aaseq1, $o_aaseq2, STDOUT);
	    
	    $aln = $factory->pairwise_alignment($seq1, $seq2);
	    $l1  = $aln->length / $o_aaseq1->length;
	    $l2  = $aln->length / $o_aaseq2->length;
	    if (($l1<$l2?$l1:$l2) > 0.7) {		
	      print "$i-$j " . $aln->percentage_identity . "% sur " . (($l1<$l2?$l1:$l2)*100) . " longueur\n";
            }
	    
	}
    }
    print "\n";
}
    
            


#$aln = $factory->pairwise_alignment($seq1, $seq2);


