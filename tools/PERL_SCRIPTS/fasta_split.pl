#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while ( my $a_ref = $fa->nextSeq ) {
  my ($name, $seq) = @{$a_ref};
  $seq = Fasta::split_large_sequence($seq, 100);
  open OUT, ">$name.fa";
  print OUT ">$name\n$seq\n";
  close OUT;
}
