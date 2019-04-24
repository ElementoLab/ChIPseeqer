#!/usr/bin/perl
use strict;



my $numtotal = 0;
my $nummapped = 0;
my $numuniq = 0;

foreach my $f (@ARGV) {
  open IN, $f or die "Cannot open file $f\n";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    next if (@a < 4);
    if ($a[2] ne "*") {
      $nummapped++;
      if ($a[11] eq "XT:A:U") {
	$numuniq++;
      }
    }
    $numtotal++;
  }
  close IN;
}


print "Total num reads = $numtotal\n";

my $r1 = sprintf("%3.1f", 100*$nummapped/$numtotal);
print "Total num mappable reads = $nummapped ($r1%)\n";

my $r2 = sprintf("%3.1f", 100*$numuniq/$numtotal);
print "Total num unambiguously mappable reads = $numuniq ($r2%)\n";
