#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];

my %H = ();

my $l = <IN>;
print $l;

while (my $l = <IN>) {
  if (!defined($H{$l})) {
    print $l;
  }
  $H{$l} = 1;
}
close IN;

