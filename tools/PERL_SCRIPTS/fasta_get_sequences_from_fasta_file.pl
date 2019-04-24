#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $h_ref = Sets::getIndex($ARGV[0]);

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  #$n = "L$n";
  if ( defined($h_ref->{$n}) ) {
    print ">$n\n$s\n\n";
  }
}
