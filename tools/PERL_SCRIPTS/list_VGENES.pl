use lib qw(/home/olly/PERL_MODULES);
use Bio::SeqIO;
use Bio::Seq;
use Sequence;



my $seqio  = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'genbank');
my $seq    = $seqio->next_seq;
my @feats = $seq->all_SeqFeatures;


foreach my $feat (@feats) {
    my $tag = $feat->primary_tag();
    if ($tag eq "V_segment") {
	print join ('', $feat->each_tag_value("gene"));
	print "\t";
	print $feat->start;
	print "\t";
	print $feat->end;
	print "\t";
	print "pseudo" if ($feat->has_tag("pseudo"));

	print "\n";
    }
}

