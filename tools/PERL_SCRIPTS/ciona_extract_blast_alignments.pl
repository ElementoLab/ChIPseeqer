use lib qw(/home/elemento/PERL_MODULES);

use Table;
use MyBlast;
use Fasta;
use Sequence;
use Repeats;
use ClustalW;
use RNAz;

use strict;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");
$mb->setDatabaseDatabase($ARGV[2]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-3");
#$mb->setMegablast(1);
#$mb->setVerbose(1);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $se = Sequence->new;
$se->setBlastDB($ARGV[1]);

my $fa = Fasta->new;


my $tmpfile1 = Sets::getTempFile("/tmp/blast.1");


my $process_tr = 0;
open LOG, ">log.txt";


foreach my $r (@$a_ref) {
  
  #print join("\t", @$r); print "\n";

  my $seq = $se->getSequenceFromBlastDB($r->[0], Sets::min($r->[3], $r->[4]), Sets::max($r->[3], $r->[4]));

  $fa->writeSeq($tmpfile1, "TT", $seq);
  $mb->setQueryDatabase($tmpfile1);
  
   my $a_ref_hsp = $mb->blastallUnique;
    
  if (scalar(@$a_ref_hsp) > 0) {
    my $h = shift @$a_ref_hsp;
#    print "$h->{QSEQ}\n$h->{DSEQ}\n";

    my @a_n = ("CI", "CS");
    my @a_s = ($h->{QSEQ}, $h->{DSEQ});
    
    my $cl = ClustalW->new;
    $cl->setSequences(\@a_n, \@a_s);
    open OUT, ">out.aln";
    print OUT $cl->getClustalWformat();
    close OUT;
    
    my $rn = RNAz->new;
    $rn->run("out.aln");

    if ($rn->isRNA()) {
      print LOG join("\t", @$r); print LOG "\t";
      print LOG "RNA\n";
      next;
      
    } else {
      print join("\t", @$r); print "\n";
    }
    
  }
}


unlink $tmpfile1;
