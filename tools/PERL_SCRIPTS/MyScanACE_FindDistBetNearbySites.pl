#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Fasta;

srand(11);

my $fa = Fasta->new;
my %L = ();
$fa->setFile($ARGV[1]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $L{$n} = length($s);
}

my %CHRT = ();
my %CHRR = ();
my %H    = ();
open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $pos = $a[1];
  my @x = ($pos, $pos+6);  
  push @{$CHRT{$a[0]}}, \@x;  

  do {
    $pos = int(rand($L{$a[0]} - 6));
  } while (defined($H{$a[0]}{$pos}));  
  $H{$a[0]}{$pos} = 1;
  my @y = ($pos, $pos+6);  
  push @{$CHRR{$a[0]}}, \@y;  

}
close IN;

my $a_ref_dt = FindDistDistances(\%CHRT);
my $a_ref_dr = FindDistDistances(\%CHRR);


for (my $i=1; $i<=25; $i++) {
  $a_ref_dt->[$i] = 0 if (!defined($a_ref_dt->[$i]));
  $a_ref_dr->[$i] = 0 if (!defined($a_ref_dr->[$i]));
  print "$i\t$a_ref_dt->[$i]\t$a_ref_dr->[$i]\n";
}



sub FindDistDistances {
  my ($h_ref_chr) = @_;
  
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
      
      $i = $j;
      next if ($d < 0);
      
      
      #print "$d\n"; 
      $D[$d] ++;
    }  
  }
  
  return \@D;
}
  

