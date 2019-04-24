#!/usr/bin/perl
#
# input : set of proteins vs set of proteins
#


#
# 1. download BLAST from http://www.ncbi.nlm.nih.gov/BLAST/download.shtml
# 2. decompress BLAST archive
# 3. modify MyBlast.pm and change the directory $self->{BLAST_DIR} = "/Users/olivier/PERL_MODULES/PROGRAMS/BLAST/bin" to the bin directory where you installed BLAST
# 4. format your input file using the BLAST utility tool formatdb
#         /path/to/blast/formatdb -i fastafile1 -p T -o T
#         /path/to/blast/formatdb -i fastafile2 -p T -o T
#

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use MyBlast;
use Fasta;
use Sets;
use Sequence;
use strict;

if (@ARGV == 0) {
  
  die "Args: fastafile1 fastafile2\n";
  
}


my $verbose = 0;
if (defined($ARGV[2])) {
  $verbose = 1;
}

my $mb = MyBlast->new;
if ($ENV{BLASTDIR} ne "") {
  $mb->setBlastDir($ENV{BLASTDIR});
}
$mb->setBlastProgram("blastp");
$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);

my $tmpfile1 = Sets::getTempFile("/tmp/blast.1");
my $tmpfile2 = Sets::getTempFile("/tmp/blast.2");

$mb->setEvalueThreshold("1e-10");
$mb->setVerbose(0);

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $qlen = undef;



#
#  get a new sequence object to retrieve BLAST seqs
#
my $s = Sequence->new;
$s->setBlastDB($ARGV[1]);
$s->useGI(0);

#
#  traverse all the proteins in file 1
#
while (my $a_ref = $fa->nextSeq) {
    
  my ($name, $seq) = @$a_ref;
  
  my $qlen = length($seq);
  
  # does not deal with small sequences
  next if ($qlen < 50);
  
  # modify name a little if we 
  if ($name =~ /^gi/) {
    my @a = split /\|/, $name;
    $name = $a[1];
  } else {
    my @a = split / /, $name;
    $name = $a[0];
  }

  if ($verbose == 1) {
    print "Starting sequence is $name\n";
  }
  
  # create a query file
  $fa->writeSeq($tmpfile1, $name, $seq);
  
  #
  # blast the sequence against the database
  #
  
  
  $mb->setDatabaseDatabase($ARGV[1]);
  $mb->setQueryDatabase($tmpfile1);
  
  my $a_ref = $mb->blastallUnique;
  
  
  next if (scalar(@$a_ref) == 0);
  
  #
  #  get the homologous protein name
  #
  my $d_id    = $mb->getUniqueHitName();
  
  if ($d_id =~ /^gi/) {
    my @a = split /\|/, $d_id;
    $d_id = $a[1];
  } 

  
  if ($verbose == 1) {
    print "Best hit is $d_id\n";
  }

  
  #
  #  get the protein sequence
  #
  if ($verbose == 1) {
    $s->setVerbose(1);
  }
  my $orth_seq = $s->getSequenceFromBlastDB($d_id, 0, 0);
  
  my $length = length($orth_seq);
  
  #
  # BLAST back
  # 

  # create a query file
  $fa->writeSeq($tmpfile2, "ORTHOLOG", $orth_seq);
  
  $mb->setBlastProgram("blastp");
  $mb->setDatabaseDatabase($ARGV[0]);
  $mb->setQueryDatabase($tmpfile2);
  
  
  my $a_ref = $mb->blastallUnique;
  
  
  my $q_id = $mb->getUniqueHitName();
  if ($q_id =~ /^gi/) {
    my @a = split /\|/, $q_id;
    $q_id = $a[1];
  } 

  if ($verbose == 1) {
    print "Best back hit is $q_id\n";
  }

  
  if ($q_id eq $name) {
    print "$name\t$d_id\n"; 
  }
  

}


unlink $tmpfile1;
unlink $tmpfile2;
