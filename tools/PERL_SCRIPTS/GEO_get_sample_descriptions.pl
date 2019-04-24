#!/usr/bin/perl

open IN, shift(@ARGV);

my @b = ("Sample_title", @ARGV);

while (my $l = <IN>) {

  chomp $l;
  foreach my $p (@b) {
    if ($l =~ /$p/) {
      my @a = split /\t/, $l, -1;
      shift @a;
      foreach my $r (@a) {
	$r =~ s/^\"//;
	$r =~ s/\"$//;	
      }
      print join("\t", @a) . "\n";
      last;
    }
  }

}
close IN;
