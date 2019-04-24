#!/usr/bin/perl

use lib "$ENV{FIREDIR}/SCRIPTS";
if ((!$ENV{FIREDIR}) || ($ENV{FIREDIR} eq '')) {
  die "Please set the FIREDIR environment variable.\n";
}

#
#  INPUT : one set of aa sequences, one list of genomes
#

use MyBlast;
use Fasta;
use Sets;
use Sequence;
use Table;

use Getopt::Long;

use strict;

my $fastafile = undef;
my $maxevalue = "1e-10";
my $nextuntil = undef;

GetOptions("fastafile=s"  =>  \$fastafile,
	   "maxevalue=s"  =>  \$maxevalue,
	   "nextuntil=s"  =>  \$nextuntil
	  );

#
#  get a new sequence object to retrieve BLAST seqs
#




my $mb = MyBlast->new;
$mb->setVerbose(1);
$mb->setBlastProgram("blastn");
$mb->setEvalueThreshold($maxevalue);
$mb->setMegablast(1);
$mb->setNbProcessors(2);
#$mb->setFilter(0);

#
#  traverse all the sequences in file
#

if (! -e "$fastafile.nin") {
  
  system("formatdb -i $fastafile -p F -o T") == 0 or die "Please format $fastafile using formatdb -i $fastafile -p F -o T\n";
  
}

my $s_tmpstore1 = Sets::getTempFile("/tmp/tmp1.seq");

my $fa = Fasta->new;
$fa->setFile($fastafile);


my $h_ref_done = undef;

#if (defined($ARGV[1])) {
#  my $ta = Table->new;
#  $ta->loadFile($ARGV[1]);
#  $h_ref_done = $ta->getIndex(0);
#}


my $skip = 1;


while (my $a_ref = $fa->nextSeq()) {
  
  my ($n, $s) = @$a_ref;

  if (defined($nextuntil) && ($n eq $nextuntil) && ($skip == 1)) {
    $skip = 0;
  }

  if ($skip == 1) {
    next;
  }

  if (defined($h_ref_done) && defined($h_ref_done->{ $n })) {
    next;
  }


  my @ss  = split /N+/, $s;
  my $max_length = 0;
  foreach my $sss (@ss) {
    if (length($sss) > $max_length) {
      $max_length = length($sss);
    }
  }
  
  if ($max_length < 28) {
    print "$n\n";
    next;
  }
    

  #my @ss  = split //, $s;
  #my @ssn = grep /N/, @ss;
  #if (scalar(@ssn) > 0.9 * scalar(@ss)) {
  #  print "$n\n";
  #  next;
  #}

  open SEQ, ">$s_tmpstore1";
  print SEQ ">$n\n$s\n\n";
  close SEQ;
  
  
  # use BLAST to align the current sequence against the reference sequence 
  $mb->setQueryDatabase($s_tmpstore1, "T");
  
  # set the database
  $mb->setDatabaseDatabase($fastafile, "T");
  
  my $a_hits = $mb->blastallMultiple;
  
  print "$n";
  
  foreach my $h (@$a_hits) {
    if ($n ne $h->{HIT_NAME}) {
      print "\t" . $h->{HIT_NAME};
      #foreach my $s ( @{ $h->{HSPS} } ) {
      #	print "  $s->{ALIGNLEN}\n";
      #      }
    }
  }

  print "\n";

}


unlink $s_tmpstore1;
