#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $a_ref = Sets::readSet($ARGV[1]);

srand;

my $n = scalar(@$a_ref);
my $a_ref_shu = Sets::shuffle_array($a_ref);

for (my $i=0; $i<$ARGV[0]; $i++) {
  print "$a_ref_shu->[$i]\n";
}

