#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

open IN, $ARGV[0] or die "Cannot open file $ARGV[0]\n";

my $l = <IN>;
chomp $l;
my @a = split /\t/, $l, -1;

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $n = shift @a;

  my @b = ();
  for (my $i=0; $i<@a; $i++) {

    if ($i == $ARGV[1]) {
      
    } else {
      push @b, $a[$i];
    }

  }

  my $avg = Sets::average(\@b);
  my $std = Sets::stddev(\@b);
  my $d   = $a[$ARGV[1]] - $avg;
  my $z   = undef;
  next if ($std == 0);
  if (abs($d) < 0.00001) {
    $z = 0;
  } else {
    $z   = $d / $std;
  }
  print "$n\t$z\n";

}
close IN;

