use lib qw(/home/elemento/PERL_MODULES);

use Table;
use MyBlast;
use Fasta;
use Sequence;
use Repeats;
use strict;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");
$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-3");
#$mb->setMegablast(1);


my $tmpfile1 = Sets::getTempFile("/tmp/blast.1");



open LOG, ">log.txt";

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {

  my ($n, $seq) = @$a_ref;
    
  $n =~ s/\-/\t/g;

  $fa->writeSeq($tmpfile1, "TT", $seq);
  $mb->setQueryDatabase($tmpfile1);
  
  my $a_ref_hits = $mb->blastallMultiple;
  
  my $cnt = 0;
  foreach my $hit (@$a_ref_hits) {
    foreach my $hsp (@{ $hit->{HSPS} }) {
      $cnt ++ if ($hsp->{EVALUE} < 1e-3);
    }
  }
  
  if ($cnt > 1) {
    
    print LOG $n; print LOG "\t";
    print LOG " > 1 hits\n";

  } else {
    
    print $n; print "\n";
    
  }
  
}

close LOG;

unlink $tmpfile1;
