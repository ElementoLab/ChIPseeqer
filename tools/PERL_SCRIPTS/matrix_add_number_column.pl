#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  $a[ $ARGV[1] ] += $ARGV[2];

  print join("\t", @a) . "\n";
  
}
close IN;

