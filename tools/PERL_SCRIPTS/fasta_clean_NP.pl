#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my @a = split /\|/, $n;
  my $r = $a[3];
  $r =~ s/\..+$//;
  next if ($r !~ /NP/);
  print ">$r\n$s\n\n";

}
