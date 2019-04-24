#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $sum = 0;
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $sum += length($s);
}

print "$sum\n";
