#!/usr/bin/perl
#
# input : set of proteins vs set of proteins
#

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use MyBlast;
use Fasta;
use Sets;
use Sequence;

my $mb = MyBlast->new;
$mb->setBlastProgram("tblastn");
$mb->setDatabaseDatabase($ARGV[1]);
my $tmpfile1 = "/tmp/blast.1";
my $tmpfile2 = "/tmp/blast.2";

$mb->setQueryDatabase($tmpfile1);
$mb->setEvalueThreshold("1e-10");
#$mb->setMismatchWeight(-);
#$mb->setVerbose(1);

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $qlen = undef;



#
#  get a new sequence object to retrieve BLAST seqs
#
my $s = Sequence->new;
$s->setBlastDB($ARGV[1]);


#
#  traverse all the proteins in file 1
#
while (my $a_ref = $fa->nextSeq) {
    
    my ($name, $seq) = @$a_ref;

    my $qlen = length($seq);

    #print "seq1=$seq\n";

    my @a = split / /, $name;
    $name = $a[0];

    # create a query file
    $fa->writeSeq($tmpfile1, $name, $seq);

    #
    # blast the sequence against the database
    #
    
    #$mb->setVerbose(1);
    #$mb->setBlastProgram("blastp");

    $mb->setFilter(0);
    $mb->setBlastProgram("blastp");
    $mb->setDatabaseDatabase($ARGV[1]);
    $mb->setQueryDatabase($tmpfile1);

    my $a_ref = $mb->blastallUnique;
    

    next if (scalar(@$a_ref) == 0);

    # put the pieces back together, from high scoring to low scoring    
    my @a_pieces = ();
    
    foreach my $r (@$a_ref) {
    	
      my $overlap = 0;
      foreach my $p (@a_pieces) {
	
	$overlap = Sets::getSequencesOverlap($r->{QFROM}, $r->{QTO}, $p->{QFROM}, $p->{QTO});
	
	last if ($overlap > 3);
      }
      
      if ($overlap <= 3) {
	push @a_pieces, $r;
      }
    }
    
    
    @a_pieces = sort { $a->{QFROM} <=> $b->{QFROM} } @a_pieces;
    
    
    #
    #  get the homologous protein name
    #
    my $d_id    = $mb->getUniqueHitName();
    
    #
    #  get the protein sequence
    #
    my $orth_seq = $s->getSequenceFromBlastDB($d_id, 0, 0);

        
    my $length = length($orth_seq);

    my $a_ref_aln = $mb->get_indentity_and_aligned_length(\@a_pieces);
    
    
    print "\%id = $a_ref_aln->[0], Laln=$a_ref_aln->[1], ql=" . $mb->getQueryLength() . " hl=$length\n";

    next unless (($a_ref_aln->[0] > 0.99) && ($a_ref_aln->[1] > 0.95 * $mb->getQueryLength()) && ($a_ref_aln->[1] > 0.95 *  $length / 3));

    #print "Total length aligned = $length, qlength = ". $mb->getQueryLength() . "\n";
    
    #print "$orth_protein\n\n";
    

    #print "seq2=$orth_protein\n";
    
    
    #
    # BLAST back
    # 
    
    
    # create a query file
    $fa->writeSeq($tmpfile2, "ORTHOLOG", $orth_seq);
    
    $mb->setFilter(1);
    $mb->setBlastProgram("blastp");
    $mb->setDatabaseDatabase($ARGV[0]);
    $mb->setQueryDatabase($tmpfile2);

    #$mb->setVerbose(1);
    
    my $a_ref = $mb->blastallUnique;


    my $q_id = $mb->getUniqueHitName();

    
    #print "hit is $q_id\n";

    if ($q_id eq $a[0]) {
	print "$a[0]\t$d_id\t$d_start\t$d_end\t$frame\n";
    }

     
    #foreach my $r (@$a_ref) {
    #print "$r->{EVALUE}\t$r->{QFROM}\t$r->{QTO}\n$r->{QSEQ}\n";
    #}
    
  }


