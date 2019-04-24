#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Fasta;
use Sets;
use MyBlast;
use strict;
use Sequence;
use Getopt::Long;

my $oligos = undef;
my $flank  = undef;
my $genome = undef;
my $match  = undef;
my $onlytop = 0;

if (@ARGV == 0) {
  die "Usage: perl map_oligos_to_genome.pl --oligos=FILE --genome=FILE --flank=1000\n";
}

GetOptions ('oligos=s'    => \$oligos,
            'genome=s'    => \$genome,
            'flank=s'     => \$flank,
            'onlytop=s'   => \$onlytop,
	    "match=s"     => \$match);


my $fa = Fasta->new;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setVerbose(0);



if (! -e "$genome.nhr") {
  $mb->setDatabaseFile($genome);
} else {
  $mb->setDatabaseDatabase($genome);
}
$mb->setMegablast(1);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-10");
$mb->setFilter(0);

#$mb->setQueryStrand(1);

    
my $se = Sequence->new;
$se->setBlastDB($genome);

#
# go thru all the sequences, align them one by one
#

mkdir "TMP" if (! -e "TMP");
my $file = Sets::getTempFile("TMP/toto");

use Table;

my $ta = Table->new;
$ta->loadFile($oligos);
my $a_ref = $ta->getArray();

print "OligoID\tseq\tchr\tmatch_start\tmatch_end\tid/aln\tstrand\te-value\tupstream ($flank bp)\tdownstream ($flank bp)\n";

foreach my $r (@$a_ref) {
  
  my $n = $r->[0];
  my $s = $r->[1];

  if (defined($match) && ($n !~ /$match/)) {
    next;
  }

  $fa->writeSeq($file, "toto", $s);
  $mb->setQueryFile($file);

  my $hsps = $mb->blastallUnique;
  my $hit  = $mb->getUniqueHitName;

  foreach my $r (@$hsps) {
    
    my $s1   = $r->{"ALIGNLEN"};
    my $e1   = $r->{"IDENTITY"};
    my $ev   = $r->{"EVALUE"};
    my $st   = $r->{"DFRAME"};    
    my $hs   = $r->{"DFROM"};
    my $he   = $r->{"DTO"};
    my $seqL = $se->getSequenceFromBlastDB($hit, $hs-$flank, $hs-1);
    my $seqR = $se->getSequenceFromBlastDB($hit, $he+1, $he+$flank);
    
    print "$n\t$s\t$hit\t$hs\t$he\t$e1/$s1\t$st\t$ev\t";

    if ($st == -1) {
      print Sets::getComplement($seqR) . "\t";
      print Sets::getComplement($seqL) . "\n";      
    } else {
      print $seqL . "\t";
      print "$seqR\n";
    }
    
    if ($onlytop == 1) {
      last;
    }
    
  }
  

}

unlink $file;
