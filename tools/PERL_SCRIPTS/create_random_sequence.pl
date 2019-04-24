#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

my @a = ("A", "C", "G", "T");

for (my $i=0; $i<$ARGV[0]; $i++) {
  my $i = rand();

  $i = int($i * 4);
  
  print $a[$i];

}

print "\n";
