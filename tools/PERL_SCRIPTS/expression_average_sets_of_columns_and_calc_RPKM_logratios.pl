#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use Sets;

use strict;



my $f  = shift @ARGV;
open IN, $f or die "cannot open $f\n"; 

my $n = scalar(@ARGV);

my $cnt = 0;
while (my $l = <IN>) {
  chomp $l;
  my @b = split /\t/, $l, -1;
  my $r = \@b;


  
  if ($cnt == 0) {

    print join("\t", @$r); #"$r->[0]";
    my $p = 1;
    for (my $i=0; $i<$n; $i++) {
      my $m = $ARGV[$i];
      my $s = "";
      my @a = ();
      for (my $j=0; $j<$m; $j++) {
	push @a, $r->[ $p ];
	$p ++;
      }
      $s = join("-", @a);
      printf("\t%s", $s);
    }
    print "\tlog2RPKMratio";
    print "\n";
    
    
  } else {
   
    my @val = ();
    my $cntval = 0;
    my $p = 1;
    for (my $i=0; $i<$n; $i++) {
      my $m = $ARGV[$i];
      my $s = 0;
      my $m_a = 0;
      for (my $j=0; $j<$m; $j++) {
	if (defined($r->[$p]) && ($r->[$p] ne "NaN") && ($r->[$p] ne "NA")) {
	  $s += $r->[ $p ];
	  $m_a ++;
	} else {
	  #die "Pb?\n";
	}
	$p ++;
      }
      if ($m_a > 0) {
	$s /= $m_a;
	push @val, sprintf("%4.3f", $s);
	$cntval ++;
      } else {
	$s = 0;
	push @val, "";
	
	#printf("\t");
      }

    }
    if ($cntval > 0) {
      my $rat = Sets::log2($val[1]+1) - Sets::log2($val[0]+1);
      print join("\t", @$r) . "\t" . join("\t", @val) . "\t$rat\n";
    }
  }
  
  $cnt ++;
}

