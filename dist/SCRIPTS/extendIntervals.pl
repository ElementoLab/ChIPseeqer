#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --intervals=FILE --ext=INT\n";
}

my $int = undef;
my $ext = 2000;

GetOptions("intervals=s"    => \$int,
           "ext=s"          => \$ext);

open IN, $int or die "Cannot open $int\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $c1 = $a[1] - $ext;
  next if ($c1 < 0);
  my $c2 = $a[2] + $ext;
  print "$a[0]\t$c1\t$c2\n";
}
close IN;

