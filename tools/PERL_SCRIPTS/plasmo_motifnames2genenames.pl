#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $a_ref = Sets::readSet($ARGV[0]);

foreach my $r (@$a_ref) {
  print Sets::Pf_motifname2genename($r) . "\n";
}
