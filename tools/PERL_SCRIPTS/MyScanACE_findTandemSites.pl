#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Fasta;

srand(11);


my %CHRT = ();
my %CHRR = ();
my %H    = ();
open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $pos = $a[1];
  my @x = ($pos, $pos+6, $a[0], $a[2], $a[3], $a[4]);  
  push @{$CHRT{$a[0]}}, \@x;  

}
close IN;

my $a_ref_dt = FindDistDistances(\%CHRT, $ARGV[1], $ARGV[2]);



sub FindDistDistances {
  my ($h_ref_chr, $tgtd1, $tgtd2) = @_;
  
  my @D = ();  
  foreach my $c (keys(%$h_ref_chr)) {  
    my @s = sort { $a->[0] <=> $b->[0] } @{$h_ref_chr->{$c}};    
    my $i = 0;
    my $n = scalar(@s);
    while ($i<$n-1) {    
      my $j = $i;
      my $d = 0;
      my $cnt = 1;
      do {
	$j ++;
	$d =  $s[$j]->[0] - $s[$i]->[1];
	#print "Try $cnt, d=$d between $s[$j]->[0]($j) and $s[$i]->[1]($i) (n=$n)\n";
	$cnt++;
      } while (($d <= 0) && ($j < $n));  
      

     
      if (($j<$n) && ($tgtd1 <= $d) && ($d <= $tgtd2)) {
	print "$d\t";
	print join("\t", @{$s[$i]}) . "\t" . join("\t", @{$s[$j]}) . "\n";
      }

      $i = $j;

      #$D[$d] ++;
    }  
  }
  
  return \@D;
}
  

