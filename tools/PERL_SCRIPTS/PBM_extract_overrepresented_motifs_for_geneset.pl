#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use strict;

use Getopt::Long;
my $geneset = undef;
my $matrix = undef;

GetOptions("geneset=s"        => \$geneset,
           "matrix=s" => \$matrix);

open IN, $matrix;

my $l = <IN>;
my @mo = split /\t/, $l;

my @idx = ();

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  

  if ($a[0] =~ /$geneset/) {
    for (my $i=0; $i<@a; $i++) {
      if ($a[$i] > 0.301) {
	push @idx, $i;
      }
    }
    last;
  }

}
close IN;

foreach my $i (@idx) {
  $mo[$i] =~ s/pwm\.txt.+$/pwm\.txt/;
  print "$mo[$i]\n";
}

