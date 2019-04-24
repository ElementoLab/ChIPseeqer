use lib qw(/home/elemento/PERL_MODULES);

use Fasta;
use Sets;
use MyBlast;
use Sequence;
use strict;



my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $se = Sequence->new;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");
#$mb->setVerbose(1);
$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-10");


my $file = Sets::getTempFile("toto");

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  
  my @a  = split /\t/, $l, -1;

  my $cid = shift @a;

  my $cnt = 0;
  my %COUNT = ();
  foreach my $eid (@a) {
    
    $eid =~ s/^r//;

    # get the cDNA with that ID
    $se->setBlastDB($ARGV[2]);
    my $seq = $se->getSequenceFromBlastDB($eid, 0, 0);
    
    next if (!$seq);

    $fa->writeSeq($file, $eid, $seq);
    $mb->setQueryFile($file);
    
    my $a_ref = $mb->blastallUnique;

    if (scalar(@$a_ref) == 0) {
	next;
    }

    my $hit_name        = $mb->getUniqueHitName;
    my $hit_length      = $mb->getUniqueHitLength;
    my $query_length    = $mb->getQueryLength;

    my $a_ref_hsps      = $mb->retain_non_overlapping_blocks($a_ref, 3);

    my $a_ref_id_le     = $mb->get_indentity_and_aligned_length($a_ref_hsps);

    my $L               = $a_ref_id_le->[1];
    my $I               = $a_ref_id_le->[0];
    
    if (($I > 0.95) && ($L > 100)) {
	#print "$cid\t$hit_name\t$L\t$I\n";
	$cnt ++;
	
	$COUNT{ $hit_name } ++;
	
	last if ($cnt == 10);
    }
  }

  if ($cnt == 0) {
    print "$cid\tNOTHING\n";
  } else {
    my $best_hit = Sets::getMaxCountFromHash(\%COUNT);
    print "$cid\t$best_hit\t$COUNT{$best_hit}\n";
  }

}

unlink $file;
