# input : set of proteins, genome
#use lib qw(/home/olly/PERL_MODULES);
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use MyBlast;
use Fasta;
use Sets;
use GeneWise;
use strict;

my $ge = GeneWise->new;

#my $verbose = 1;

my $mb = MyBlast->new;
$mb->setBlastProgram("tblastn");
$mb->setDatabaseDatabase($ARGV[1]);
my $tmpfile1 = Sets::getTempFile("blast.1");
my $tmpfile2 = Sets::getTempFile("blast.2");
my $tmpfile3 = Sets::getTempFile("blast.3");

$mb->setQueryDatabase($tmpfile1);
$mb->setEvalueThreshold("1e-2");
$mb->setNbProcessors(2);
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
if (defined($ARGV[3])) {
  $h_ref_todo = Sets::getIndex($ARGV[3]);
}

open OUT1, ">$ARGV[2].seq2.aa";
open OUT2, ">$ARGV[2].seq2.nt";    

while (my $a_ref = $fa->nextSeq) {
    
    my ($name, $seq) = @$a_ref;

    next if (defined($h_ref_todo) && !defined($h_ref_todo->{ $name }));

    #if ($verbose == 1) {
    #  print "$name\n";
    #}

    my $qlen = length($seq);

    my @a = split / /, $name;

    # create a query file
    $fa->writeSeq($tmpfile1, $name, $seq);

    # do the blast
    #$mb->setVerbose(1);
   
    $mb->setBlastProgram("tblastn");
    $mb->setDatabaseDatabase($ARGV[1]);
    $mb->setQueryDatabase($tmpfile1);

    my $a_ref = $mb->blastallUnique;
    
    next if (scalar(@$a_ref) == 0);

    # put the pieces back together, from high scoring to low scoring
        
    # init the first cluster
    my @CLUSTERS    = ();
    my $cnt_cluster = 0;
    
    $CLUSTERS[ $cnt_cluster ]->{ DFROM } = $a_ref->[0]->{ DFROM };
    $CLUSTERS[ $cnt_cluster ]->{ DTO   } = $a_ref->[0]->{ DTO   };
    $CLUSTERS[ $cnt_cluster ]->{ DFRAME} = $a_ref->[0]->{ DFRAME} / abs($a_ref->[0]->{ DFRAME});
    $a_ref->[0]->{ TAKEN } = 1;

    #print "seed:\n";
    #print "  $a_ref->[0]->{EVALUE}\n";    
    #print "  $a_ref->[0]->{DFROM }\n";    
    #print "  $a_ref->[0]->{DTO   }\n";    
    #print "  $a_ref->[0]->{DFRAME}\n";

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



    my $re_fr = $CLUSTERS[ $cnt_cluster ]->{ DFROM } - 200;
    $re_fr = 1 if ($re_fr < 1);
    my $re_to = $CLUSTERS[ $cnt_cluster ]->{ DTO   } + 200;
    $re_to = $l_id if ($re_to > $l_id);
    # 
    #  get the big fragment
    #
    my $d_id    = $mb->getUniqueHitName();


    my $region = $o_seq->getSequenceFromBlastDB($d_id, ($re_fr<1?1:$re_fr), $re_to);

    if ($CLUSTERS[ $cnt_cluster ]->{ DFRAME } < 0) {
      $region = Sets::getComplement($region);
    }


    #print "region extracted = $re_fr -> $re_to\n";
    
    #print "$region\n";
    $fa->writeSeq($tmpfile2, "ORTHOLOG", $region);

    #
    #  use genewise between query and region
    #
    
    my $a_res = $ge->run($tmpfile1, $tmpfile2);

    

    #print join("\t", @$a_res); print "\n";

    my $l_offset  = $a_res->[0] - 1;
    my $r_offset  = length($region) - $a_res->[1] ;

    #print "$l_offset\t$r_offset\n";
    
    my $d_start = undef;
    my $d_end   = undef;
    my $d_frame = $CLUSTERS[ $cnt_cluster ]->{ DFRAME };
    if ($CLUSTERS[ $cnt_cluster ]->{ DFRAME } > 0) {
      $d_start = $re_fr + $l_offset;
      $d_end   = $re_to - $r_offset;
    } else {
      $d_start = $re_fr + $r_offset;
      $d_end   = $re_to - $l_offset;
    }
    
    if ($d_start < 1) {
      $d_start = 1;
    }
    
    #
    # get the exons
    #
    my $a_ref_exons = $ge->getExons();
    my $genomic_seq = "";
    foreach my $e (@$a_ref_exons) {
      $genomic_seq .= substr($region, $e->[0]-1, $e->[1]-$e->[0]+1);
    }


    #
    # BLAST back
    # 

    # create a query file
    
    $fa->writeSeq($tmpfile3, "ORTHOLOG", $a_res->[3]);
    $mb->setBlastProgram("blastp");
    $mb->setDatabaseDatabase($ARGV[0]);
    $mb->setQueryDatabase($tmpfile3);
    my $a_ref = $mb->blastallUnique;
    my $q_id = $mb->getUniqueHitName();
    if ($q_id eq $a[0]) {
      print "$a[0]\t$d_id\t$d_start\t$d_end\t$d_frame\n";

      print OUT1 ">$a[0]\n$a_res->[3]\n\n";
      print OUT2 ">$a[0]\n$genomic_seq\n\n";

    } 
      
      
}

close OUT1;
close OUT2;

unlink $tmpfile1;
unlink $tmpfile2;
unlink $tmpfile3;
