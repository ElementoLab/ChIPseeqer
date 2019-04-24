use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;
use Sets;
use MyBlast;
use strict;

my $adalen = 45;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setVerbose(1);
#$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("0.01");
$mb->setFilter(0);
$mb->setGapOpening(2);
$mb->setGapExtension(1);
$mb->setq(-1);
#$mb->setWordLength(100);
$mb->setQueryDatabase($ARGV[0]);

#my $a_ref_reads = $mb->blastallMultiple_Unique(100);

open IN, $ARGV[0]; my @lines = <IN>; close IN; 
my $txt = join("", @lines);

my $a_ref_reads = $mb->_analyzeMultiple_UniqueXML($txt, 100);
foreach my $r (@$a_ref_reads) {


  my $numhsps = @{$r->{HSPS}};
  my $hsp     = $r->{HSPS};

 

  
  my $i1 = Sets::min($hsp->[0]->{QFROM}, $hsp->[0]->{QTO});
  my $i2 = Sets::max($hsp->[0]->{QFROM}, $hsp->[0]->{QTO});

  
  my $dup = 0;
  for (my $i=1; $i<$numhsps; $i++) {
    my $j1 = Sets::min($hsp->[$i]->{QFROM}, $hsp->[$i]->{QTO});
    my $j2 = Sets::max($hsp->[$i]->{QFROM}, $hsp->[$i]->{QTO});
    
    my $d1 = Sets::getSequencesOverlap($i1, $i2,
				       $j1, $j2);     
    
    if (($d1 > 0) && ($hsp->[0]->{EVALUE} == $hsp->[$i]->{EVALUE})) {
      $dup = 1;
    }
  }  
  if ($dup == 1) {
    #print "\tduplicate";
    next;
  }
  
  # else best match is ok
  # print hsp0
  my $txt = "$r->{QUERY}->{NAME}\t$r->{QUERY}->{LENGTH}\t$r->{HIT}->{NAME}";
  if ($ARGV[1] ne "") {
    $txt = "$ARGV[1]\t$txt";
  }
  print "$txt";
  print "\t$hsp->[0]->{EVALUE}";
  print "\t$hsp->[0]->{QFROM}";
  print "\t$hsp->[0]->{QTO}";
  print "\t$hsp->[0]->{QFRAME}";
  print "\t$hsp->[0]->{DFROM}";
  print "\t$hsp->[0]->{DTO}";
  print "\t$hsp->[0]->{DFRAME}";
  print "\t-\t-\t-";
  my $dahsp = $mb->getHSPdata($hsp->[0]->{QSEQ}, $hsp->[0]->{DSEQ});
  print "\t" . join("\t", @$dahsp);
  print "\n";

  # check dist on genome
  for (my $i=1; $i<$numhsps; $i++) {
    
    my $j1 = Sets::min($hsp->[$i]->{QFROM}, $hsp->[$i]->{QTO});
    my $j2 = Sets::max($hsp->[$i]->{QFROM}, $hsp->[$i]->{QTO});

    my $d1 = Sets::getSequencesOverlap($i1, $i2,
				       $j1, $j2);     
    
    my $d2 = Sets::getSequencesOverlap($hsp->[0]->{DFROM}, $hsp->[0]->{DTO},
				      $hsp->[$i]->{DFROM}, $hsp->[$i]->{DTO});     
   
    my $lab = undef;
    if (($d2 > 0) && ($d1 < 0) && (abs($d1) >= 0.5*$adalen) && (abs($d1) <= 1.5*$adalen)) {
      $lab = "SMRTbell-separated";
      # distances are within 20% of each other 
    } elsif (($d1 < 0) && ($d2 < 0) && ($hsp->[0]->{QFRAME} == $hsp->[$i]->{QFRAME}) && ($hsp->[0]->{DFRAME} == $hsp->[$i]->{DFRAME}) && ((abs($d1)/abs($d2)) > 0.75) && ((abs($d1)/abs($d2)) < 1.25)) {
      $lab = "Compatible";
    } else {
      next;
    }

    print "$txt";
    print "\t$hsp->[$i]->{EVALUE}";
    print "\t$hsp->[$i]->{QFROM}";
    print "\t$hsp->[$i]->{QTO}";
    print "\t$hsp->[$i]->{QFRAME}";
    print "\t$hsp->[$i]->{DFROM}";
    print "\t$hsp->[$i]->{DTO}";
    print "\t$hsp->[$i]->{DFRAME}";
    print "\t$d1\t$d2\t$lab";
    my $dahsp = $mb->getHSPdata($hsp->[$i]->{QSEQ}, $hsp->[$i]->{DSEQ});
    print "\t" . join("\t", @$dahsp);
    print "\n";
 
  }

}
