#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

if (@ARGV == 0) {
  die "Args: matrix col1 col2 (ratios=col1/col2)[ threshold]\n";
}
open IN, $ARGV[0] or die "Cannot open file";

my $i1 = $ARGV[1];
my $i2 = $ARGV[2];

my $l = <IN>; chomp $l;

my @a = split /\t/, $l;

if ($ARGV[3] eq "") {
  print "TS\t$a[$i1]/$a[$i2]\n";
}
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  
  my $lr = Sets::log2( ($a[$i1]+1) / ($a[$i2]+1) );
  
  if ($ARGV[3] eq "") {
    print "$a[0]\t$lr\n";
  } elsif (($ARGV[3] > 0) && ($lr >= Sets::log2($ARGV[3]))) {
    print "$a[0]\n";
  } elsif (($ARGV[3] < 0) && ($lr <= Sets::log2(-1/$ARGV[3]))) {
    print "$a[0]\n";
  }

}
close IN;

