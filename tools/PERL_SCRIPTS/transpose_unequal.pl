#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my $n = @$a_ref;

my @M = ();
my $maxm = 0;

for (my $i=0; $i<$n; $i++) {

  if (@{$a_ref->[$i]} > $maxm) {
    $maxm = @{$a_ref->[$i]};
  }

  for (my $j=0; $j<@{$a_ref->[$i]}; $j++) {
    $M[$j][$i] = $a_ref->[$i]->[$j];
  }

}


for (my $i=0; $i<$maxm; $i++) {
  
  my @a = ();
  for (my $j=0; $j<2; $j++) {
    push @a, $M[$i][$j];
  }
  print join("\t", @a) . "\n";
}
