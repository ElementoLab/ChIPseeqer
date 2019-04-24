#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

open IN, $ARGV[0];

my $l = <IN>;

print "ID_REF\tMAX\n";

while (my $l = <IN>) {  
  chomp $l;  
  my @a = split /\t/, $l, -1;
  my $p = shift @a;
  my $m = Sets::maxInArray(\@a); 
  print "$p\t$m\n";
  
}

close IN;


