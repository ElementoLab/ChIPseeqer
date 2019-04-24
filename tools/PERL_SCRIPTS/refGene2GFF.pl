#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

open IN, $ARGV[0];

my $genome = $ARGV[1]; #"mm9_refGene";

if (!defined($ARGV[1]) || ($ARGV[1] eq "")) {
  $genome = "hg18_refGene";
}

my %TSS = ();
my %NM  = ();

while (my $l = <IN>) {

  chomp $l;
  my @a = split /\t/, $l, -1;

  my $n      = $a[1];
  #next if ($n ne "NM_001164233");
  $NM{$n} ++ ;

  if ($NM{$n} > 1) {
    my $nn = $NM{$n}-1;
    $n = "$n\_dup$nn";
  }

  my $c      = $a[2];
  my $fr     = $a[3];
  my $ts_st  = $a[4]+1;
  my $ts_en  = $a[5];
  my $cds_st = $a[6]+1;
  my $cds_en = $a[7];
  my $e_st   = $a[9];
  my $e_en   = $a[10];
  my $e_fr   = $a[15];
  my $g      = $a[12];

  my $tss    = undef;
  if ($fr eq "+") {
    $tss = $ts_st;
  } else {
    $tss = $ts_en;
  }

  # add TSS if not in there yet
  if (!Sets::in_array(@{$TSS{$g}}, $tss)) {
    push @{$TSS{$g}}, $tss;
  }
  
  my $idx = -1;
  my $cnt = 1;
  foreach my $t (@{$TSS{$g}}) {
    if ($t == $tss) {
      $idx = $cnt;
    }
    $cnt++;
  }
  if ($idx == -1) {
    die "Problem\n";
  }

  # desc
  my $desc = "gene_id \"$g\"; transcript_id \"$n\"; g_id \"$g\"; tss_id \"$g.$idx\";";

  my $tes    = undef;
  if ($fr eq "+") {
    $tes = $ts_en;
  } else {
    $tes = $ts_st;
  }

  $e_st =~ s/\,$//;
  $e_en =~ s/\,$//;
  $e_fr =~ s/\,$//;

  my @a_e_st = split /\,/, $e_st;
  my @a_e_en = split /\,/, $e_en;
  my @a_e_fr = split /\,/, $e_fr;

  if ($fr eq "+") {
    
  }

  for (my $i=0; $i<@a_e_st; $i++) {
    my $ex_st = $a_e_st[$i]+1;
    my $ex_en = $a_e_en[$i];

    my $ex_fr = $a_e_fr[$i];
    if ($ex_fr == 1) {
      $ex_fr = 2;
    } elsif ($ex_fr == 2) {
      $ex_fr = 1;
    }

    if (($n =~ /^NM/) && Sets::sequencesOverlap($ex_st, $ex_en, $cds_st, $cds_en)) {

      my $ex_cds_st = $ex_st;
      if (Sets::sequencesOverlap($ex_st, $ex_en, $cds_st, $cds_st)) {  # if exon overlaps with CDS start
	$ex_cds_st = $cds_st;
	
	if ($fr eq "+") { # start
	  #my $start_end = $ex_cds_st + 2;
	  #print "$c\t$genome\tstart_codon\t$ex_cds_st\t$start_end\t0.000000\t$fr\t.\t$desc\n";
	} else { #stop
	  #my $start_end = $ex_cds_st + 2;
	  #print "$c\t$genome\tstop_codon\t$ex_cds_st\t$start_end\t0.000000\t$fr\t.\t$desc\n";
	}
	
      }
      my $ex_cds_en = $ex_en;
      if (Sets::sequencesOverlap($ex_st, $ex_en, $cds_en, $cds_en)) {  # if exon overlaps with CDS start
	$ex_cds_en = $cds_en;
	
	if ($fr eq "+") {   # stop codon
	  #$ex_cds_en  -= 3;                # walk back two
	  #my $stop_en  = $ex_cds_en + 1;   
	  #my $stop_beg = $ex_cds_en + 3;
	  #print "$c\t$genome\tstop_codon\t$stop_en\t$stop_beg\t0.000000\t$fr\t.\t$desc\n";
	} else { # start codon	  
	  #$ex_cds_en  -= 3;                # walk back two
	  #my $stop_en  = $ex_cds_en + 1;   
	  #my $stop_beg = $ex_cds_en + 3;
	  #print "$c\t$genome\tstart_codon\t$stop_en\t$stop_beg\t0.000000\t$fr\t.\t$desc\n";

	}
	

      }
      print "$c\t$genome\tCDS\t$ex_cds_st\t$ex_cds_en\t0.000000\t$fr\t$ex_fr\t$desc\n";
    }
    

  
    
    print "$c\t$genome\texon\t$ex_st\t$ex_en\t0.000000\t$fr\t.\t$desc\n";

    
    #chr1    mm9_refGene     start_codon     134212807       134212809       0.000000        +       .       gene_id "Nuak2"; transcript_id "NM_028778"; g_id "Nuak2"; tss_id "Nuak2.1";
    #chr1    mm9_refGene     CDS     134212807       134213049       0.000000        +       0       gene_id "Nuak2"; transcript_id "NM_028778"; g_id "Nuak2"; tss_id "Nuak2.1";
    #chr1    mm9_refGene     exon    134212702       134213049       0.000000        +       .       gene_id "Nuak2"; transcript_id "NM_028778"; g_id "Nuak2"; tss_id "Nuak2.1";


  }

}

close IN;

