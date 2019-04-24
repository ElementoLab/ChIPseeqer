#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

open IN, $ARGV[0];
print "Motif 1\n";

my $n = undef;
while (my $l = <IN>) {
  chomp $l;
  $n = length($l);
  
  print "$l\n";
  

}
close IN;

print "*" x $n;
print "\n";



