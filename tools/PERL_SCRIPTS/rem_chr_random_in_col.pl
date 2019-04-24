#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $v = $a[$ARGV[1]];

  if (($v !~ /random/) && ($v !~ /hap/)) {
    print "$l\n";
  } 

}
close IN;

