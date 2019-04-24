use lib qw(/home/elemento/PERL_MODULES);

use Table;
use MyBlast;
use Fasta;
use Sequence;
use Repeats;
use strict;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");
$mb->setDatabaseDatabase($ARGV[2]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-3");
$mb->setMegablast(1);
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
  
  

  my $seq = $se->getSequenceFromBlastDB($r->[0], Sets::min($r->[3], $r->[4]), Sets::max($r->[3], $r->[4]));

  if ($process_tr == 1) {
    # first run TRF
    my $tr = Repeats->new;
    $tr->setSeq($seq);
    $tr->setProgram('TRF');
    $tr->process();
    my $a_ref_tandems = $tr->getResults;
    $tr->dispose;
    
    my $seq_masked = Sets::maskExons($seq, $a_ref_tandems, 'X');

    my @ax  = split //, $seq_masked;
    my @axr = grep /X/, @ax;
    
    if (scalar(@axr) / scalar(@ax) > 0.5) {
      print LOG join("\t", @$r); print LOG "\t";
      print LOG "tandem repeat\n";
      next;
    }
  }

  $fa->writeSeq($tmpfile1, "TT", $seq);
  $mb->setQueryDatabase($tmpfile1);
  
   my $a_ref_hsp = $mb->blastallUnique;
    
  if (scalar(@$a_ref_hsp) > 0) {
    
    my $txt = $mb->getUniqueHitName();

    my $a_ref_il = $mb->get_indentity_and_aligned_length($a_ref_hsp);
    

    print LOG join("\t", @$r); print LOG "\t";
    print LOG "match RNA $txt, i=$a_ref_il->[0], l=$a_ref_il->[1]\n";

  } else {
    
    print join("\t", @$r); print "\n";
    #print "no RNA match\n";
    
  }

  #<STDIN>;
  
}

close LOG;

unlink $tmpfile1;
