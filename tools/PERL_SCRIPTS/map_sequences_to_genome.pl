BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Fasta;
use Sets;
use MyBlast;
use strict;



my $fa = Fasta->new;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setVerbose(0);
$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-20");
$mb->setMegablast(1);
$mb->setFilter(0);
$mb->setWordLength(100);


#$mb->setQueryStrand(1);

    
    
#
# go thru all the sequences, align them one by one
#

my $file = Sets::getTempFile("toto");

$fa->setFile($ARGV[0]);


print "Seq\tchr\tmatch_start\tmatch_end\tstrand\tqlen\tid/aln\te-value\n";
while (my $a_ref_mirna = $fa->nextSeq()) {

  my $n = $a_ref_mirna->[0]; $n =~ s/\ .+$//;
  my $s = $a_ref_mirna->[1];

  my $l = length($s);

  $fa->writeSeq($file, $n, $s);
  $mb->setQueryFile($file);

  my $hsps = $mb->blastallUnique;
  
  if (@$hsps > 0) {

    # get best hsps
    my $r = shift @$hsps;
    
    my $ch = $mb->getUniqueHitName();
    my $s1 = $r->{"ALIGNLEN"};
    my $e1 = $r->{"IDENTITY"};
    my $ev = $r->{"EVALUE"};
    my $st = $r->{"DFRAME"};
    
    my $hs = $r->{"DFROM"};
    my $he = $r->{"DTO"};
    
    # num identities in alignment / length aligned region
    my $fr = $e1 / $s1;  
    
    if (($s1 >= 0.99*$l) && ($fr >= 0.99)) {    
      print "$n\t$ch\t$hs\t$he\t$st\t$l\t$e1/$s1\t$ev\n";
    } else {
      print "FAILED\t$n\t$ch\t$hs\t$he\t$st\t$l\t$e1/$s1\t$ev\n";
    }

  }
  

}

unlink $file;
