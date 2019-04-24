#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
#my $l = <IN>; chomp $l;
#my @a = split /\t/
my $cnt = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $n = shift @a;
  print "$n";
  for (my $i=0; $i<@a-1; $i++) {
    for (my $j=$i+1; $j<@a; $j++) {
      if ($cnt > 0) {
	my $r = log($a[$j]+1) - log($a[$i]+1);
	print sprintf("\t%3.2f", $r);
      } else {
	print "\t$a[$j]/$a[$i]";
      }
    }
  }
  print "\n";
  $cnt ++;
}
close IN;

