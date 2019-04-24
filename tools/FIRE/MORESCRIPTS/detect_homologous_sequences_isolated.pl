

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
use strict;

my $df = DataFiles->new;

#
#  get a new sequence object to retrieve BLAST seqs
#


my $mb = MyBlast->new;
$mb->setVerbose(0);
$mb->setBlastProgram("blastn");
$mb->setEvalueThreshold("1e-10");
$mb->setMegablast(1);
$mb->setNbProcessors(2);
$mb->setFilter(0);
$mb->setQueryStrand(1);  # top strand
#
#  traverse all the sequences in file
#



my $s_tmpstore1 = Sets::getTempFile("/tmp/tmp1.seq");

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my $h_ref_done = undef;
if (defined($ARGV[1])) {
  my $ta = Table->new;
  $ta->loadFile($ARGV[1]);
  $h_ref_done = $ta->getIndex(0);
}


while (my $a_ref = $fa->nextSeq()) {
  
  my ($n, $s) = @$a_ref;

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
  $mb->setDatabaseDatabase($ARGV[0], "T");
  
  my $a_hits = $mb->blastallMultiple;

  
  print "$n";
  

  foreach my $h (@$a_hits) {
    if ($n ne $h->{HIT_NAME}) {

      my $really = 0;
      foreach my $s ( @{ $h->{HSPS} } ) {
	#print " * l=$s->{ALIGNLEN} *";
      	if ($s->{ALIGNLEN} >= 200) {
	  $really = 1; last;
	}
      }

      if ($really == 1) {
	print "\t" . $h->{HIT_NAME};
      }

    }
  }

  print "\n";

}


unlink $s_tmpstore1;
