use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;
use Sets;
use MyBlast;
use strict;
use  ClustalW;

my $aligner = "clustalw";  # dialign

my $adalen = 45;


# load reads
my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
my %READS = ();
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $READS{$n} = $s;
}
$fa->dispose();


# load miRNAs
$fa = Fasta->new;
$fa->setFile($ARGV[1]);
my %REFS = ();
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $REFS{$n} = $s;
}
$fa->dispose();




my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setVerbose(1);
$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("0.01");
#$mb->setFilter(0);
$mb->setGapOpening(2);
$mb->setGapExtension(1);
#$mb->setq(-1);
#$mb->setWordLength(100);
$mb->setQueryDatabase($ARGV[0]);

my $a_ref_reads = $mb->blastallMultiple_Unique(100);

#open IN, $ARGV[0]; my @lines = <IN>; close IN; 
#my $txt = join("", @lines);

#my $a_ref_reads = $mb->_analyzeMultiple_UniqueXML($txt, 100);

my %MAP = ();

foreach my $r (@$a_ref_reads) {


  my $numhsps = @{$r->{HSPS}};
  my $hsp     = $r->{HSPS};

 

  
  #my $i1 = Sets::min($hsp->[0]->{QFROM}, $hsp->[0]->{QTO});
  #my $i2 = Sets::max($hsp->[0]->{QFROM}, $hsp->[0]->{QTO});

  
  #my $dup = 0;
  #for (my $i=1; $i<$numhsps; $i++) {
  #  my $j1 = Sets::min($hsp->[$i]->{QFROM}, $hsp->[$i]->{QTO});
  #  my $j2 = Sets::max($hsp->[$i]->{QFROM}, $hsp->[$i]->{QTO});
  #  
  #  my $d1 = Sets::getSequencesOverlap($i1, $i2,
  #$j1, $j2);     
  #
  #  if (($d1 > 0) && ($hsp->[0]->{EVALUE} == $hsp->[$i]->{EVALUE})) {
  #    $dup = 1;
  #  }
  #}
  #if ($dup == 1) {
  #print "\tduplicate";
  #next;
  #}
  
  # else best match is ok
  # print hsp0
  my $txt = "$r->{QUERY}->{NAME}\t$r->{QUERY}->{LENGTH}\t$r->{HIT}->{NAME}";
  #if ($ARGV[1] ne "") {
  #  $txt = "$ARGV[1]\t$txt";
  #}

  my $qlenfrac = $hsp->[0]->{ALIGNLEN} / $r->{QUERY}->{LENGTH};
  my $hlenfrac = $hsp->[0]->{ALIGNLEN} / $r->{HIT}->{LENGTH};


  print "$txt";
  print "\t$hsp->[0]->{EVALUE}";
  print "\t$hsp->[0]->{ALIGNLEN}";
  print "\t$hsp->[0]->{QFROM}";
  print "\t$hsp->[0]->{QTO}";
  print "\t$hsp->[0]->{QFRAME}";
  print "\t$hsp->[0]->{DFROM}";
  print "\t$hsp->[0]->{DTO}";
  print "\t$hsp->[0]->{DFRAME}";
  print "\t-\t-\t-";
  my $dahsp = $mb->getHSPdata($hsp->[0]->{QSEQ}, $hsp->[0]->{DSEQ});
  print "\t" . join("\t", @$dahsp);

  my $skip = 0;
  if (($qlenfrac < 0.50) || ($hlenfrac < 0.50)) {
    #next;
    print "\tTOO SHORT";
    $skip = 1;
  }

  print "\n";

  next if ($skip == 1);

  # store
  my $seq = $READS{$r->{QUERY}->{NAME}};
  if ($hsp->[0]->{DFRAME} < 0) {
    $seq = Sets::getComplement($seq);
  }
  my @a_seq = ($r->{QUERY}->{NAME}, $seq);
  push @{$MAP{$r->{HIT}->{NAME}}}, \@a_seq;

  # check dist on genome
  #for (my $i=1; $i<$numhsps; $i++) {
    
  #  my $j1 = Sets::min($hsp->[$i]->{QFROM}, $hsp->[$i]->{QTO});
  #  my $j2 = Sets::max($hsp->[$i]->{QFROM}, $hsp->[$i]->{QTO});
  #
  #  my $d1 = Sets::getSequencesOverlap($i1, $i2,
  #$j1, $j2);     
  #  
  #  my $d2 = Sets::getSequencesOverlap($hsp->[0]->{DFROM}, $hsp->[0]->{DTO},
  #$hsp->[$i]->{DFROM}, $hsp->[$i]->{DTO});     
   
  #  my $lab = undef;
  #  if (($d2 > 0) && ($d1 < 0) && (abs($d1) >= 0.5*$adalen) && (abs($d1) <= 1.5*$adalen)) {
  #    $lab = "SMRTbell-separated";
  #    # distances are within 20% of each other 
  #  } elsif (($d1 < 0) && ($d2 < 0) && ($hsp->[0]->{QFRAME} == $hsp->[$i]->{QFRAME}) && ($hsp->[0]->{DFRAME} == $hsp->[$i]->{DFRAME}) && ((abs($d1)/abs($d2)) > 0.75) && ((abs($d1)/abs($d2)) < 1.25)) {
  #    $lab = "Compatible";
  #  } else {
  #    next;
  #  }
  #
  #  print "$txt";
  #  print "\t$hsp->[$i]->{EVALUE}";
  #  print "\t$hsp->[$i]->{QFROM}";
  #  print "\t$hsp->[$i]->{QTO}";
  #  print "\t$hsp->[$i]->{QFRAME}";
  #  print "\t$hsp->[$i]->{DFROM}";
  #  print "\t$hsp->[$i]->{DTO}";
  #  print "\t$hsp->[$i]->{DFRAME}";
  #  print "\t$d1\t$d2\t$lab";
  #  my $dahsp = $mb->getHSPdata($hsp->[$i]->{QSEQ}, $hsp->[$i]->{DSEQ});
  #  print "\t" . join("\t", @$dahsp);
  #  print "\n";
 
  #}

}



foreach my $m (keys(%MAP)) {
  
  print "$m\n";
  my $numseqs = 0;
  my $txt = "";
  foreach my $s (@{$MAP{$m}}) {
    print "$s->[1]\n";
    $txt .= ">$s->[0]\n$s->[1]\n";
    $numseqs ++;
  }
  $txt .= ">$m\n$REFS{$m}\n";
  print "\n";
  
  
  
  print "Writing merged file for $m.\n";
  my $ffa = "$ARGV[0].$m";
  open OUT, ">$ffa";
  print OUT "$txt";
  close OUT;



  if ($numseqs <= 200) {
    
    print "Running dialign ... ";
    
    #  -thr 5
    my $todo = "$ENV{HOME}/PERL_MODULES/PROGRAMS/DIALIGN/dialign2-2 -fa -n $ffa";
    print "$todo\n";
    system($todo);
    
    # transform
    open OUT, ">$ffa.aln";
    my $fa = Fasta->new;
    $fa->setFile("$ffa.fa");  
    while (my $a_ref = $fa->nextSeq()) {
      my ($n, $s) = @$a_ref;
      print OUT "$n\t$s\n";
    }
    $fa->dispose();
    close OUT;

    unlink "$ffa.ali";
    unlink "$ffa.fa";

  } else {
    
    my $todo = "clustalw -INFILE=$ffa -TYPE=DNA -OUTFILE=$ffa.clu -OUTORDER=INPUT > /dev/null ";
    print "$todo\n";
    system($todo) == 0 or die "Cannot exec $todo\n";

    my $cl = ClustalW->new;
    my $numchar = length("GJP4OI302IZGP0-F7-R6      ");
    $cl->setNumCharName($numchar);        
    $cl->setFile("$ffa.clu");
    my $a_ref_aln = $cl->getSeqsWithNames();
    open OUT, ">$ffa.aln";
    foreach my $s (@$a_ref_aln) {
      print  OUT "$s->[0]\t$s->[1]\n";
    }
    close OUT;

    
  }

  
  print "Done.\n"; 

}
