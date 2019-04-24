#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;
use Sets;
use MyBlast;
use strict;



my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-3");
#$mb->setQueryStrand(1);
$mb->setVerbose(1);
$mb->setGapOpening(2);
$mb->setGapExtension(1);
$mb->setq(-1);
    

$mb->setQueryFile($ARGV[0]);
    
my $a_ref = $mb->blastallMultiple;


foreach my $hit (@$a_ref) { 
  
  my $hsps = $hit->{"HSPS"};
  
  @$hsps = sort { $a->{DFROM} <=> $b->{DFROM} } @$hsps;
  
  foreach my $r (@$hsps) {
    
    my $s1 = $r->{"ALIGNLEN"};
    my $e1 = $r->{"IDENTITY"};
    my $ev = $r->{"EVALUE"};
    my $st = $r->{"DFRAME"};
    my $po = $r->{"DFROM"};
    print "$hit->{HIT_NAME}\t$e1\t$ev\t$st\t$po\n";
  }
  print "\n";
}

