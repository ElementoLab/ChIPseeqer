#
# input : set of proteins vs set of proteins
#

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use lib "$home/usr/lib/perl5/site_perl/5.6.1";
use lib "$home/usr/lib/perl5/site_perl/5.8.3";

use MyBlast;
use Fasta;
use Sets;
use Sequence;
use Table;
use Bio::SearchIO;
use DataFiles;
use strict;

my $df = DataFiles->new;

#
#  get a new sequence object to retrieve BLAST seqs
#
my $s = Sequence->new;
$s->setBlastDB($ARGV[1]);


my $mb = MyBlast->new;
$mb->setBlastProgram("blastp");
$mb->setEvalueThreshold("1e-10");
$mb->setDatabaseDatabase($ARGV[1]);
$mb->setStoreSequences(1);


#
#  traverse all the proteins in file 1
#
my $a_ref_genes    = Sets::readSet($ARGV[0]);

my $nbgenes     = scalar(@$a_ref_genes);
my $s_tmpstore1 = Sets::getTempFile("/tmp/tmp1.seq");

for (my $i=0; $i<$nbgenes; $i++) {


    #
    # get the protein from DB and save it into a temp file
    #
    my $seq1 = $s->getSequenceFromBlastDB($a_ref_genes->[$i], 0, 0);
    open SEQ, ">$s_tmpstore1";
    print SEQ ">$a_ref_genes->[$i]\n$seq1\n\n";
    close SEQ;
    
    print ">$a_ref_genes->[$i]\n$seq1\n";
    
    
    # use BLAST to align the current sequence against the reference sequence 
    $mb->setQueryFile($s_tmpstore1, "T");

    my $a_hits = $mb->blastallMultiple;

    

    my $query_length = $mb->getQueryLength();

    foreach my $myhit (@$a_hits) {

	my $nbhsps     =  scalar(@{ $myhit->{"Hit_hsps"} });
	my $hit_name   =  $myhit->{"Hit_id"};
	
	#
	#  if the current hit is an isoform, skip it
	#
	if ($hit_name eq $a_ref_genes->[$i]) {
	    #print " isoform\n";
	    next;
	}
	
	
	my $hit_length = $myhit->{"Hit_len"};
	
	
	#
	#  order all HSPs, remove overlapping ones
	#
	@{ $myhit->{"Hit_hsps"} } = sort { $a->{"Hsp_query-from"} <=> $b->{"Hsp_query-from"} } @{ $myhit->{"Hit_hsps"} };
	my @a_hsps = ();
	foreach my $h (@{ $myhit->{"Hit_hsps"} }) {
	    my $overlap = 0;
	    foreach my $p (@a_hsps) {
		$overlap = Sets::getSequencesOverlap($h->{"Hsp_query-from"}, $h->{"Hsp_query-to"}, 
						     $p->{"Hsp_query-from"}, $p->{"Hsp_query-to"});
		last if ($overlap > 3);
	    }
	    
	    if ($overlap <= 3) {
		push @a_hsps, $h;
	    }
	}
	@a_hsps = sort { $a->{"Hsp_query-from"} <=> $b->{"Hsp_query-from"} } @a_hsps;
	

	#
	#  calculate the matching length, and the identity
	#
	my $l1 = 0;
	my $l2 = 0;
	my $gapped_seq1 = "";
	my $gapped_seq2 = "";

	foreach my $h (@a_hsps) {
	    my $s1 =  $h->{"Hsp_qseq"};  	    
	    $gapped_seq1 .= $s1;
	    
	    $s1 =~ s/\-//g;
	    $l1 += length($s1);
	    
	    my $s2 =  $h->{"Hsp_hseq"};  
	    $gapped_seq2 .= $s2;

	    $s2 =~ s/\-//g;
	    $l2 += length($s2);
	}
	
	#  compute the identity between the natching segments
	my @a           = split //, $gapped_seq1;
        my @b           = split //, $gapped_seq2;
        my $len         = scalar(@a);
        my $cnt_matches = 0;
        my $cnt_gaps1   = 0;
        my $cnt_gaps2   = 0;
        for (my $j=0; $j<$len; $j++) {
            $cnt_matches ++ if ($a[$j] eq $b[$j]);
            $cnt_gaps1   ++ if ($a[$j] eq '-');
            $cnt_gaps2   ++ if ($b[$j] eq '-');
        }

	#print "SEQ1=$gapped_seq1\n";
	#print "SEQ2=$gapped_seq2\n";


        #
        #   calculate the real fraction of identical residues
        #
	
	#print "L=" . length($gapped_seq1) . "\n";
	#print "$cnt_gaps1 - $cnt_gaps2\n";

        my $frac_identical = ($cnt_matches / ( length($gapped_seq1) - $cnt_gaps1 - $cnt_gaps2 ));

        #
        #   conserved fraction in each sequence
        #
        my $frac_len_seq1    = ( length($gapped_seq1) - $cnt_gaps1 ) / $query_length;
        my $frac_len_seq2    = ( length($gapped_seq2) - $cnt_gaps2 ) / $hit_length;

        #
        #  are these two genes paralogous ?
        #
        my $paralog = 0;
        my $L = length($gapped_seq1) - ( $cnt_gaps1 + $cnt_gaps2 ) / 2;
        my $I = $frac_identical;
        my $IMAX = undef;
        my $LMIN = 0.5 * Sets::max($query_length, $hit_length);  # get 50% of the length of the largest protein

        if ($L >= 150) {
            $IMAX = 0.30;
        } else {

            $IMAX = 0.01 * 6 + 4.8 * $L**(- 0.32 * (1 + exp(- $L / 1000.0)));
        }

        if (($I >= $IMAX) && ($L > $LMIN)) {
            $paralog = 1;
        }


        if ($paralog == 1) {
	    print $a_ref_genes->[$i] . "\t" . $myhit->{"Hit_id"} . "\n";
	    #print "ORTHOLOGS !\n";
	}
	
    }
}


unlink $s_tmpstore1;

