#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  if ($s !~ /\./) {
    print ">$n\n$s\n";
  }

}
