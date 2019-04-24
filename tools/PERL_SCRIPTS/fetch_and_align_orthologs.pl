#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use strict;

#1. read all human sequences NP
#2. read datasets, 
#3. for each dataset, read seqs, for each orths, add seq to human seq 
# fill a table gene x species => seq

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my %MAT = ();
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $MAT{$n}{"Human"} = $s; 
}
$fa->dispose();

open IN, "datasets.txt" or die "nods\n";
while (my $l = <IN>) {
  chomp $l;
  my ($sp,$da) = split /[\t\ ]+/, $l, -1;
  my $ff = Sets::filename($da);
  $ff =~ s/\.gz//;
#print "opening $ff\n";
  # read fasta
  my $fi = Fasta->new;
  $fi->setFile($ff);  
  my %SEQ = ();
  while (my $a_ref = $fi->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $n =~ s/\ .+$//;
    $SEQ{$n}= $s; 
     #print "SEQ{$n}=$s\n";
  }
  $fi->dispose();

  # read orthologs
  open IN2, "$sp.txt" or die "Cannot open $sp.txt\n";
  while (my $l = <IN2>) {
    chomp $l;
    my ($np,$or) = split /\t/, $l, -1;
    $MAT{$np}{$sp} = $SEQ{$or};
  }
  close IN2;

  
  
}
close IN;

my $dir = "ALN";
mkdir $dir if (! -e $dir);
foreach my $np (keys(%MAT)) {
  
  open OUT, ">$dir/$np.fa";
  foreach my $sp (keys(%{$MAT{$np}})) {
    print OUT ">$sp\n$MAT{$np}{$sp}\n";
  }
  close OUT;

}
