#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";


open IN, $ARGV[0] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[2] > $a[1]) {
    print "$a[0]\t$a[1]\t$a[2]\n";
  } else {
    print "$a[0]\t$a[2]\t$a[1]\n";
  }
}
close IN;

