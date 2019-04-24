#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my $l = $ARGV[2] - $ARGV[1] + 1;

  my $ss = substr($s, $ARGV[1]-1, $l);
  
  print ">$n\n$ss\n";
}
