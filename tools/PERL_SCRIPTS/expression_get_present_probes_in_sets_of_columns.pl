

use lib qw(/home/elemento/PERL_MODULES);
use strict;



my $f  = shift @ARGV;
open IN, $f or die "cannot open $f\n"; 

my $n = scalar(@ARGV);

my $cnt = 0;
while (my $l = <IN>) {
  chomp $l;
  my @b = split /\t/, $l, -1;
  my $r = \@b;

  print "$r->[0]";
  
  if ($cnt == 0) {
    my $p = 1;
    for (my $i=0; $i<$n; $i++) {
      my $m = $ARGV[$i];
      my $s = "";
      my @a = ();
      for (my $j=0; $j<$m; $j++) {
	push @a, $r->[ $p ];
	$p ++;
      }
      $s = join("/", @a);
      printf("\t%s", $s);
    }
    print "\n";
    
    
  } else {
   
    my $p = 1;
    my $P_row = 0;
    for (my $i=0; $i<$n; $i++) {
      my $m = $ARGV[$i];
      my $s = 0;
      my $m_a = 0;
      my %H_P_LOCAL = ();
      for (my $j=0; $j<$m; $j++) {
	$H_P_LOCAL{$r->[$p]} ++;
	$p ++;
      }
      
      if ($H_P_LOCAL{'P'} > $m/2) {
	$P_row ++;
      }
    }
    
    print "\t$P_row\n";
    
  }
  
  $cnt ++;
}

