#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


open IN, $ARGV[0] or die "Cannot open file";

my $l = <IN>; chomp $l;

my @a = split /\t/, $l;

print "TS\tEXP\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  
  my $lr = Sets::log2( ($a[3]+1) / ($a[2]+1) );
  
  if ($ARGV[1] eq "") {
    print "$a[0]\t$lr\n";
  } elsif (($ARGV[1] > 0) && ($lr >= Sets::log2($ARGV[1]))) {
    print "$a[0]\n";
  } elsif (($ARGV[1] < 0) && ($lr <= Sets::log2(-1/$ARGV[1]))) {
    print "$a[0]\n";
  }

}
close IN;

