# input : set of proteins, genome
#use lib qw(/home/olly/PERL_MODULES);
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use MyBlast;
use Fasta;
use Sets;
use Sim4;
use strict;

my $si = Sim4->new;

my $verbose = 0;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");
$mb->setDatabaseDatabase($ARGV[1]);
my $tmpfile1 = Sets::getTempFile("blast.1");
my $tmpfile2 = Sets::getTempFile("blast.2");
my $tmpfile3 = Sets::getTempFile("blast.3");

$mb->setQueryDatabase($tmpfile1);
$mb->setEvalueThreshold("1e-2");
$mb->setNbProcessors(2);
$mb->setFilter(1);
$mb->setExitOnCrash(0);


#$mb->setMismatchWeight(-);
#$mb->setVerbose(1);

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $qlen = undef;

my $o_seq = Sequence->new;
$o_seq->setBlastDB($ARGV[1]);

#
#  load a list of sequence for parallelization purpose ...
#
my $h_ref_todo = undef;
if (defined($ARGV[2])) {
  $h_ref_todo = Sets::getIndex($ARGV[2]);
}

#open OUT1, ">$ARGV[2].seq2.aa";
#open OUT2, ">$ARGV[2].seq2.nt";    


#my $start = 0;
while (my $a_ref = $fa->nextSeq) {
    
    my ($name, $seq) = @$a_ref;


    #print "$name ... ";
    
    #if ($name =~ /XPF2n0796/) {
    #  $start = 1;
    #}
    
    #if ($start == 0) {
    #print " ignored.\n";
    #next;
    #} else {
    #print " continue.\n";
    #}

    # next if (defined($h_ref_todo) && !defined($h_ref_todo->{ $name }));

    #if ($verbose == 1) {
    #  print "$name\n";
    #}

    my $qlen = length($seq);

    my @a = split / /, $name;

    # create a query file
    $fa->writeSeq($tmpfile1, $name, $seq);

    # do the blast
    #$mb->setVerbose(1);
   
    $mb->setBlastProgram("blastn");
    $mb->setDatabaseDatabase($ARGV[1]);
    $mb->setQueryDatabase($tmpfile1);

    my $a_ref = $mb->blastallUnique;
    
    next if (($mb->crashed()) || (scalar(@$a_ref) == 0));

    # put the pieces back together, from high scoring to low scoring
        
    # init the first cluster
    my @CLUSTERS    = ();
    my $cnt_cluster = 0;
    
    $CLUSTERS[ $cnt_cluster ]->{ DFROM } = $a_ref->[0]->{ DFROM };
    $CLUSTERS[ $cnt_cluster ]->{ DTO   } = $a_ref->[0]->{ DTO   };
    $CLUSTERS[ $cnt_cluster ]->{ DFRAME} = $a_ref->[0]->{ DFRAME} / abs($a_ref->[0]->{ DFRAME});
    $a_ref->[0]->{ TAKEN } = 1;

    # eliminate the highest e-value
    shift @$a_ref;   

    # sort by increasing dist on chr
    @$a_ref = sort { $a->{DFROM} <=> $b->{DFROM} } @$a_ref;

    #print $CLUSTERS[ $cnt_cluster ]->{DFRAME} . "$ev\t\n";
    foreach my $r (@$a_ref) {
      $r->{DFRAME} = $r->{DFRAME} / abs( $r->{ DFRAME } );
    }
    
    # get the cluster containing the highest e-value, basically ...
    my $changed = 1;
    while ($changed == 1) {

      # traverse all fragments
      $changed = 0;
	    
      foreach my $r (@$a_ref) {

	#print "add fragment $r->{DFROM} -> $r->{DTO}, $r->{DFRAME}";
	#print "  $r->{EVALUE}\n";    


	# if a fragment is not part of a cluster, consider it
	if (!defined( $r->{ TAKEN } ) && ($r->{DFRAME} == $CLUSTERS[ $cnt_cluster ]->{ DFRAME })) {
	  
	  # calculate the distance between the two HSPs
	  my $overlap = Sets::getSequencesOverlap($r->{DFROM},
						  $r->{DTO}  ,
						  $CLUSTERS[ $cnt_cluster ]->{ DTO   },
						  $CLUSTERS[ $cnt_cluster ]->{ DFROM });
	  # if distance < 100000, add to cluster
	  if ($overlap > -100000) {
	    $CLUSTERS[ $cnt_cluster ]->{ DTO   } = Sets::maxInArray( [ $r->{DFROM},
								       $r->{DTO}  ,
								       $CLUSTERS[ $cnt_cluster ]->{ DTO   },
								       $CLUSTERS[ $cnt_cluster ]->{ DFROM } ] );
	    $CLUSTERS[ $cnt_cluster ]->{ DFROM } = Sets::minInArray( [ $r->{DFROM},
								       $r->{DTO}  ,
								       $CLUSTERS[ $cnt_cluster ]->{ DTO   },
								       $CLUSTERS[ $cnt_cluster ]->{ DFROM } ] );


	    #print "sure\n";
	    $r->{TAKEN} = 1;
	    $changed = 1;
	  } else {
	    #print "no, too far\n";
	  }
	} #end if
      } # end foreach
    } #end while

   # print $CLUSTERS[ $cnt_cluster ]->{ DFROM   } . "\t" . $CLUSTERS[ $cnt_cluster ]->{ DTO } . "\n";


    my $l_id    = $mb->getUniqueHitLength();


    my $ext     = 1000;

    my $re_fr = $CLUSTERS[ $cnt_cluster ]->{ DFROM } - $ext;
    $re_fr = 1 if ($re_fr < 1);

    # real offset
    my $real_ext = $CLUSTERS[ $cnt_cluster ]->{ DFROM } - $re_fr;

    my $re_to = $CLUSTERS[ $cnt_cluster ]->{ DTO   } + $ext;
    $re_to = $l_id if ($re_to > $l_id);

    #
    #  get the big fragment
    #
    my $d_id    = $mb->getUniqueHitName();
    my $region  = $o_seq->getSequenceFromBlastDB($d_id, ($re_fr<1?1:$re_fr), $re_to);

    my $strand  = $CLUSTERS[ $cnt_cluster ]->{ DFRAME };

    #if ($CLUSTERS[ $cnt_cluster ]->{ DFRAME } < 0) {
    #  $region = Sets::getComplement($region);
    #}


    if ($verbose == 1) {
      print "region extracted for sim4 = [ $re_fr -> $re_to ]\n";
    }

    #print "$region\n";
    $fa->writeSeq($tmpfile2, "ORTHOLOG", $region);

    #
    #  use sim4 between query and region for a fine exon stricture match
    #

    if ($verbose == 1) {
      print "Running sim4.\n";
    }
    my $h_ref_data = $si->run($tmpfile1, $tmpfile2);

    foreach my $e (@{$h_ref_data->{RAW}}) {
      #print join("\t", @$e) . "\n";
    }
    

    my $len_matched = undef;
    if (($h_ref_data->{MINC} == 1) && ($h_ref_data->{MAXC} == $qlen)) {
      #print "cDNA matches genome across its entire length.\n";
      $len_matched = 100;
    } else {
      
      $len_matched = ($h_ref_data->{MAXC} - $h_ref_data->{MINC} + 1) / $qlen;
      $len_matched = int(0.5 + 100*$len_matched);

    }

    #print "Aligning $name.\n";

    my @nn = split /\|/, $name;

    $name = pop @nn;
    $name =~ s/\(/\t/;
    $name =~ s/\)//;

    if ($verbose == 1) {
      print "cDNA matches across $len_matched% length, with $h_ref_data->{AVGSIM}% average similarity.\n", ;
    }

    my $gen_st = $re_fr + $h_ref_data->{MING} - 1;
    my $gen_en = $re_fr + $h_ref_data->{MAXG} - 1;
    
    #print "Genomic region matches = [ $gen_st, $gen_en ]\n";
    my $num_exons = scalar(@{$h_ref_data->{RAW}});

    print "$name\t$d_id\t$gen_st\t$gen_en\t$strand\t$num_exons\t$len_matched\t$h_ref_data->{AVGSIM}\n";

}

#close OUT1;
#close OUT2;

unlink $tmpfile1;
unlink $tmpfile2;
unlink $tmpfile3;
