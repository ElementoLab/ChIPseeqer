#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use strict;
use Sets;

my $a_ref = Sets::get_kmers_matching_re($ARGV[0]);

foreach my $r (@$a_ref) {
  print "$r\n";
}
