#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use strict;


open IN, $ARGV[0];
my %H = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $H{"$a[0]\t$a[1]\t$a[2]"} = $a[3];
}
close IN;

open IN, $ARGV[1];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $peak = "$a[0]\t$a[1]\t$a[2]";
  if (defined($H{$peak})) {
    my $rat = $a[3] - $H{$peak};     
    print "$peak\t$rat\n";
  } else {
    print "$peak not defined\n";
  }

}
close IN;



