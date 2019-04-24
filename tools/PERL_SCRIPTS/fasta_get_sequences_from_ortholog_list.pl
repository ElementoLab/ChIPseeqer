#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES /home/elemento/PERL_MODULES);

use Sequence;
use Sets;
use Table;
use Getopt::Long;
use strict;

if (scalar(@ARGV) == 0) {
    die "Usage : perl fasta_get_sequences_from_ortholog_list.pl --orthologs=FILE --fasta1=FILE --fasta2=FILE --outfasta1=FILE --outfasta2=FILE\n";
}

my $orthologs = undef;
my $fasta1    = undef;
my $fasta2    = undef;
my $outfasta1 = undef;
my $outfasta2 = undef;


GetOptions ('orthologs=s'  => \$orthologs,
	    'fasta1=s'     => \$fasta1,
	    'fasta2=s'     => \$fasta2,
	    'outfasta1=s'  => \$outfasta1,
	    'outfasta2=s'  => \$outfasta2,
	    
	   );

open OUT1, ">$outfasta1" or die "Cannot open $outfasta1.\n";
open OUT2, ">$outfasta2" or die "Cannot open $outfasta2.\n";

system("formatdb -i $fasta1 -p F -o T") if (! -e "$fasta1.nhr");
system("formatdb -i $fasta2 -p F -o T") if (! -e "$fasta2.nhr");

open INCH, $orthologs or die "Cannot open $orthologs.\n";

if (! -e "/home/olly/PERL_MODULES/PROGRAMS/BLAST/bin/fastacmd") {
  die "fastacmd not found. Do a locate fastacmd, then update this script with the correct directory.\n";
}

while (my $l = <INCH>) {
  chomp $l;
  my @a = split /\t/, $l;
  
  my $l1 = $a[0];
  my $l2 = (defined($a[1])?$a[1]:$a[0]);	

  print "Processing $l1 / $l2 .. ";
  
  my $o_seq = Sequence->new;
  #$o_seq->setVerbose(1);
  $o_seq->setBlastPath("/home/olly/PERL_MODULES/PROGRAMS/BLAST/bin");
  
  $o_seq->setBlastDB($fasta1);    
  my $s1 = $o_seq->getSequenceFromBlastDB("$l1", 0, 0);
  
  $o_seq->setBlastDB($fasta2);    
  my $s2 = $o_seq->getSequenceFromBlastDB("$l2", 0, 0);

  if (defined($s1) && ($s1 ne "") && defined($s2) && ($s2 ne "")) {
    print " both defined.\n";
    print OUT1 ">$l1\n$s1\n\n";
    print OUT2 ">$l2\n$s2\n\n";
  } else {
    print " not defined.\n";
  }
  
}  
    
    
close INCH;
close OUT1;
close OUT2;
