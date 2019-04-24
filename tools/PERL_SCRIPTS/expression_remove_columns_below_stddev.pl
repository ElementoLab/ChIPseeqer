#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;
my @V = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $n = shift @a;
  for (my $i=0; $i<@a; $i++) {
    push @{$V[$i]}, $a[$i];
  }
 
}
close IN;

my @SD = ();
for (my $i=0; $i<@V; $i++) {
  $SD[$i] = Sets::stddev($V[$i]);
}

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $n = shift @a;
  print $n;
  for (my $i=0; $i<@a; $i++) {
    if ($SD[$i] > $ARGV[1]) {
      print "\t$a[$i]";
    }
  }
  print "\n";
}
close IN;

