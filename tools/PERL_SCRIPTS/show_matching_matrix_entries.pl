#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;


my $i=0; 

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $r = \@a;
  for (my $j=0; $j<@$r; $j++) {
    if ($r->[$j] =~ /$ARGV[1]/) {
      print "$i\t$j\t$r->[$j]\n";
    }
  }
  
  $i ++;
}

close IN;
