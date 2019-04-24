#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my @a = split //, $s;
  
  my $cntstop = 0;
  for (my $i=0; $i<@a-1; $i++) {
    if ($a[$i] eq '*') {
      $cntstop ++;
    }
  }
  
  if ($cntstop == 0) {
    print ">$n\n$s\n\n";
  }

}
