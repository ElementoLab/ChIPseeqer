#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Graph;

my $gr = Graph->new;

$gr->loadGraph($ARGV[0]);

my $a_ref = $gr->getCenteredGraph($ARGV[1], $ARGV[2]);


foreach my $r (@$a_ref) {
  print join("\t", @$r) . "\n";
}
