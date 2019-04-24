#
#  INPUT : one set of aa sequences, one list of genomes
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
use DataFiles;
use File::Basename;
use strict;


my $df = DataFiles->new;

#
#  get a new sequence object to retrieve BLAST seqs
#


my $mb = MyBlast->new;
#$mb->setVerbose(1);
$mb->setBlastProgram("blastp");
$mb->setEvalueThreshold("1e-10");
$mb->setNbProcessors(2);

my $a_ref_genomes    = Sets::readSet($ARGV[1]);

#
#  if a list of genes is defined .. use it.
#
my $h_ref_genes = undef;
if (defined($ARGV[2])) {
  $h_ref_genes = Sets::getIndex($ARGV[2]);
}


my $nbgenomes     = scalar(@$a_ref_genomes);
my $tmpfile1 = Sets::getTempFile("/tmp/tmp1.seq");
my $tmpfile2 = Sets::getTempFile("/tmp/tmp2.seq");

#
#  traverse all the proteins in file 1

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $se = Sequence->new;



while (my $a_ref = $fa->nextSeq()) {
    
    my ($n, $s) = @$a_ref;

    #next if (length($s) < 50);

    #
    # get gene name 
    #
    my ($gn) = $n =~ /(.+?)$/;

    if (defined($h_ref_genes) && (!defined($h_ref_genes->{ $gn }))) {
      next;
    }

    print "$gn";
    
    #
    # get the protein from DB and save it into a temp file
    #
    open SEQ, ">$tmpfile1";
    print SEQ ">$n\n$s\n\n";
    close SEQ;
    
    #
    #  go thru all the genomes
    #

    for (my $i=0; $i<$nbgenomes; $i++) {

	next if (! -e "$a_ref_genomes->[$i].psd");
	
	#
	#  exclude own genome
	#
	next if (basename($ARGV[0]) eq basename($a_ref_genomes->[$i]));

	
	# extract the nameof the genome
	my ($g) = $a_ref_genomes->[$i] =~ /GENOMES\/(.+?)\/genome/;
	
	# set the query and database
	$mb->setQueryDatabase($tmpfile1, "T");
	$mb->setDatabaseDatabase($a_ref_genomes->[$i], "T");

	# run the blast
	my $a_hits = $mb->blastallMultiple;
	my $query_length = $mb->getQueryLength();
	
	my $cnt = 0;
        if (scalar(@$a_hits) > 0) { 

	  #print "got " . scalar(@$a_hits) . " hits\n";
	  
	  my $cnt_processed_hits = 0;
	  foreach my $hit (@$a_hits) {

	    $cnt_processed_hits++;
	    
	    #last if ($cnt_processed_hits == 6);
	    
	    my $a_ref_hsps = $hit->{HSPS};
	    my $hit_length = $hit->{HIT_LENGTH};
	    
	    #
	    #  order all HSPs, remove overlapping ones
	    #
	    @$a_ref_hsps = sort { $a->{QFROM} <=> $b->{QFROM} } @$a_ref_hsps;
	  	    
	    $a_ref_hsps = $mb->retain_non_overlapping_blocks($a_ref_hsps, 3);
	    
	    #
	    #  calculate the matching length, and the identity
	    #
	    my $a_ref_il = $mb->get_indentity_and_aligned_length($a_ref_hsps);
	    
	    my ($I, $L) = @$a_ref_il;
	    
	    my $LMIN = (2.0 / 3.0) * Sets::max($query_length, $hit_length);  # get 50% of the length of the largest protein
	    	    
	    if ($L >= $LMIN) {

	      #
	      #  do the reciprocal BLAST here ..
	      #
	      
	      

	      
	      # get the protein sequence
	      my $orth_protein = "";
	      #foreach my $r (@$a_ref_hsps) {
	      #	$orth_protein .= $r->{DSEQ};
	      #      }
	      #      $orth_protein =~ s/\-//g;
	      
	      # alternative
	      my $hit_name = $hit->{HIT_NAME}; $hit_name =~ s/\ +//g;
	      #print "got a hit, $hit_name ..\n";
	      $se->setBlastDB($a_ref_genomes->[$i]);
	      #$se->setVerbose(1);
	      $orth_protein = $se->getSequenceFromBlastDB($hit_name, 0, 0);
	      
	      #print "got $orth_protein\n";

	      
	      
	      
	      $fa->writeSeq($tmpfile2, "ORTHOLOG", $orth_protein);
	      
	      $mb->setBlastProgram("blastp");
	      $mb->setDatabaseDatabase($ARGV[0]);
	      $mb->setQueryDatabase($tmpfile2);
	      
	      #
	      #  two cases: blast against own genome or NOT
	      #
	      if (basename($ARGV[0]) ne basename($a_ref_genomes->[$i])) {
		
		my $a_ref = $mb->blastallUnique;
		
		my $q_id = $mb->getUniqueHitName();
		$q_id =~ s/\ +//g;
		
		#print "reci BLAST returned $q_id\n";
		
		if ($q_id eq $gn) {
		  $cnt++;
		} else {
		  last;
		}
		
	      } # end:dealing with other genome
	      #
	      #  own genonme
	      #
	      else {
		
		# look at the second best-reciprocal hit (or the first)
		my $a_hits = $mb->blastallMultiple;
		
		#
		# if we get at least one hit
	        #
		if (scalar(@$a_hits) > 0) {
		  
		  # get the first, do nothing for now
		  my $hit1      = shift @$a_hits; 
		  my $hit1_name = $hit1->{HIT_NAME}; $hit1_name =~ s/\ +//g;
		  
		  # get the second
		  if (scalar(@$a_hits) > 0) {
		    
		    my $hit2      = shift @$a_hits; 
		    my $hit2_name = $hit2->{HIT_NAME}; $hit2_name =~ s/\ +//g;
		    
		    # if the second hit is reciprocal, inc the counter
		    if ($hit2_name eq $gn) {
		      $cnt ++;

		      # else maybe the first hit is reciprocal; weird but also increment the counter
		    } elsif ($hit1_name eq $gn) {
		      $cnt ++;
		    }
		    
		    # we have only one hit, but it is reciprocal, so inc the counter
		  } elsif ($hit1_name eq $gn) {
		    $cnt++;
		  }
		  
		}
	      } # end: dealing with own genome
	      
	    } 
	    # no appropriate math -> finish
	    else {
	      last;
	    } 
	  }

	}

	#$mb->dispose();   # free memory

	print " $cnt";

	
	
    }

    print "\n";
}
     
	
unlink $tmpfile1;
unlink $tmpfile2;
