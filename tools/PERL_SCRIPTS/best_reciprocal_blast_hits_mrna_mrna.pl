#!/usr/bin/perl
#
# input : set of mrna vs mrna
#

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use MyBlast;
use Fasta;
use Sets;
use Sequence;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");
$mb->setDatabaseDatabase($ARGV[1]);
my $tmpfile = "/tmp/blast.1";
$mb->setQueryDatabase($tmpfile);
$mb->setEvalueThreshold("1e-20");
$mb->setNbProcessors(2);
#$mb->setMegablast(1);

#$mb->setVerbose(1);

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $qlen = undef;



#
#  get a new sequence object to retrieve BLAST seqs
#
my $s = Sequence->new;
$s->setBlastDB($ARGV[1]);


my $h_ref_done = undef;
if (defined($ARGV[2])) {
  my $ta = Table->new;
  $ta->loadFile($ARGV[2]);
  $h_ref_done = $ta->getIndex(0);
}

#
#  traverse all the proteins in file 1
#
while (my $a_ref = $fa->nextSeq) {
    
    my ($name, $seq) = @$a_ref;

    if (defined($h_ref_done) && defined($h_ref_done->{ $name })) {
      next;
    }

    my $qlen = length($seq);

    my @a = split / /, $name;
    $name = $a[0];

    # create a query file
    $fa->writeSeq($tmpfile, $name, $seq);

    #
    # blast the sequence against the database
    #
    
    #$mb->setVerbose(1);
    #$mb->setBlastProgram("blastp");
    
    $mb->setDatabaseDatabase($ARGV[1]);
    $mb->setQueryDatabase($tmpfile);

    my $a_ref = $mb->blastallUnique;
    if (scalar(@$a_ref) == 0) {
      print "$a[0]\tnot found\n";
      next;
    }

    #
    #  get the homologous protein name
    #
    my $d_id    = $mb->getUniqueHitName();
    
    #
    #  get the protein sequence
    #
    my $orth_seq = $s->getSequenceFromBlastDB($d_id, 0, 0);

    #
    # BLAST back
    # 

    # create a query file
    $fa->writeSeq($tmpfile, "ORTHOLOG", $orth_seq);

    $mb->setDatabaseDatabase($ARGV[0]);
    $mb->setQueryDatabase($tmpfile);
    
    my $a_ref = $mb->blastallUnique;

    my $q_id = $mb->getUniqueHitName();


    if ($q_id eq $a[0]) {
      print "$a[0]\t$d_id\n";
    } else {
      print "$a[0]\tnot reciprocal\n";
    }

}


