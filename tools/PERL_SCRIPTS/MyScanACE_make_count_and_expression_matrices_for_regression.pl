#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use Getopt::Long;
use strict;

my $motifdir = undef;
my $expdata  = undef;
my $genelist = "/Users/olivier/PROGRAMS/PLASMODIUM/GENOMES/annotations/pf_annotation_6.0.txt.genenames";
my $mult     = 0;
my $weight   = 0;
my $suffix   = undef;
my $mingenes = 100;

if (@ARGV == 0) {
  print "Args: --motifdir=DIR --expdata=FILE\n";
  exit;
}

GetOptions("motifdir=s" => \$motifdir,
	   "genelist=s" => \$genelist,
	   "weight=s"   => \$weight,
	   "mult=s"     => \$mult,
	   "suffix=s"   => \$suffix,
           "expdata=s"  => \$expdata);


# get gene list from sequences
my $a_ref_seqgenes = Sets::readSet ($genelist);
my $h_ref_seqgenes = Sets::getIndex($genelist);

#
# read motif matches
#
opendir(DIR, $motifdir) || die "can't opendir $motifdir: $!";
my @files = grep { /\.matches\.txt$/ } readdir(DIR);
closedir DIR;


my %COUNTS = ();
foreach my $f (@files) {

  my $maxscore = 0;
  if ($weight == 1) {
    # first pass to determine max score
    open IN, "$motifdir/$f" or die "cannot open $f\n";
    while (my $l = <IN>) {
      chomp $l;
      my @a = split /\t/, $l, -1;
      if ($a[3] > $maxscore) {
	$maxscore = $a[3];
      }
    }
    print STDERR "max score for $f is $maxscore\n";
    close IN;
  }

  open IN, "$motifdir/$f" or die "cannot open $f\n";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    # motif x gene matrix
    my $add = 1;
    if ($weight == 1) {
      $add = $a[3] / $maxscore;
      #print STDERR "$l, Adding score = $add\n";
    }

    $COUNTS{$f}{$a[0]} += $add;

  }
  close IN;
}

#
# read exp data
#


# first pass to eliminate motifs with not enough
my %MATCHCOUNT = ();
my @MOTIFS     = keys(%COUNTS);
open IN, $expdata or die "cannot open $expdata\n";
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if (defined($h_ref_seqgenes->{$a[0]})) {
    
    # loop over mo
    foreach my $m (@MOTIFS) {      
      if (defined($COUNTS{$m}{$a[0]}) && ($COUNTS{$m}{$a[0]} > 0)) {
	$MATCHCOUNT{$m} ++;
      } else {
	$COUNTS{$m}{$a[0]} = 0;
      }
      
    }
	  
  }

}
close IN;

my @NEWMOTIFS = ();
foreach my $m (keys(%MATCHCOUNT)) {
  if ($MATCHCOUNT{$m} >= $mingenes) {
    push @NEWMOTIFS, $m;
  } else {
    print STDERR "Skipping $m\n";
  }
}
@MOTIFS = @NEWMOTIFS;

print STDERR scalar(@NEWMOTIFS) . " with num gene > $mingenes\n";


#
# OUTPUT
# 
my $oute = "$expdata.expcount.exp";
my $outc = "$expdata.expcount.mot";
if (defined($suffix)) {
  $oute .= $suffix;
  $outc .= $suffix;
}

open OUTE, ">$oute";
open OUTC, ">$outc";

#
# HEADER
#
my @MOTH   = @MOTIFS;
map s/\_pwm.*$//, @MOTH;

print OUTE $l;
print OUTC "MOTIFS\t" . join("\t", @MOTH);
if ($mult == 1) {
  my $n = @MOTH;
  for (my $i=0; $i<$n-1; $i++) {
    for (my $j=$i+1; $j<$n; $j++) {
      print OUTC "\t$MOTH[$i]:$MOTH[$j]";
    }
  }
}
print OUTC "\n";


#
# COUNTS
#
open IN, $expdata or die "cannot open $expdata\n";
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if (defined($h_ref_seqgenes->{$a[0]})) {
    print OUTE "$l\n";
    print OUTC "$a[0]";
    foreach my $m (@MOTIFS) {
      #if (!defined($COUNTS{$m}{$a[0]})) {
      #$COUNTS{$m}{$a[0]} = 0;
      #}
      print OUTC "\t" . $COUNTS{$m}{$a[0]};
    }
    # add columns for mult
    if ($mult == 1) {
      my $n = @MOTIFS;
      for (my $i=0; $i<$n-1; $i++) {
	for (my $j=$i+1; $j<$n; $j++) {
	  print OUTC "\t" . $COUNTS{$MOTIFS[$i]}{$a[0]} * $COUNTS{$MOTIFS[$j]}{$a[0]};
	}
      }
    }
    print OUTC "\n";
  }

}
close IN;


close OUTC;
close OUTE;

print "$oute created\n";
print "$outc created\n";
