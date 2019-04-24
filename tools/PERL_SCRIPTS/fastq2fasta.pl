#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;



open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";

my $outfile = "$ARGV[0].fa";
my $outqual = "$outfile.qual";
open OUT1, ">$outfile" or die "Cannot open $outfile\n";
open OUT2, ">$outqual" or die "Cannot open $outqual\n";

my $pth = 1;
if ($ARGV[1] ne "") {
  $pth = $ARGV[1];
}
my $minlen = 0;
if ($ARGV[2] ne "") {
  $minlen = $ARGV[2];
}
my $prefix = "";
if ($ARGV[3] ne "") {
  $prefix = $ARGV[3];
}

while (my $l = <IN>) {
  chomp $l;
  #my @a = split /\t/, $l, -1;
  
  if ($l =~ /^\@/) {    
    my ($num) = $l =~ /\/(\d+)$/;
    my $l2 = <IN>; chomp $l2;
    my $l3 = <IN>;
    my $l4 = <IN>; chomp $l4;
    my $len = length($l2);
    my @b = split //, $l4;

    if (@b == 0) {
      #die ("Odd: $l\n");
      next;
    }
    
    my $sum = 0;
    foreach my $s (@b) {
      my $qs = ord($s); 
      $qs -= 33;
      my $qsp = exp( $qs / (-10 * log(10)));

      $sum += $qsp;
      #print "\n";
    }
    $sum /= @b;

    my $p = $sum;

    #$sum -= 33; 
    #my $p = exp( $sum / (-10 * log(10)));    
    #printf "Q=%3.2f, p=%3.2f\n", $sum, $p;
    #<STDIN>;
    
    if (($p <= $pth) && ($len >= $minlen)) {
      print OUT1 ">read$num\n$l2\n";
      #$p = sprintf("%3.2f", $p);
      print OUT2 "$prefix\tread$num\t$len\t$p\n";
    }
  }

}
close IN;
close OUT1;
close OUT2;
