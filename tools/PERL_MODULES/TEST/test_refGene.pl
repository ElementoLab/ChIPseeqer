#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use refGene;


my $re = refGene->new;

my $a_ref = $re->FindGenesThatOverlapWith("chr17", 7512444, 7531588);

foreach my $r (@$a_ref) {
  print join("\t", @$r) . "\n";
}
