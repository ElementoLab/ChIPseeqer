#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  #$n = "L$n";
  if ($n =~ /$ARGV[1]/) { 
    $s = uc($s);
    print ">$n\n$s\n\n";
  }
}
