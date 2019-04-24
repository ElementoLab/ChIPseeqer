#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use strict;



my $f  = shift @ARGV;
open IN, $f or die "cannot open $f\n"; 

my $n = scalar(@ARGV);

my $cnt = 0;
my $l   = <IN>;
print $l;
while (my $l = <IN>) {
  chomp $l;
  my @b = split /\t/, $l, -1;
  my $r = \@b;

  my @val = ();
  my $cntval = 0;
  my $p = 1;
  for (my $i=0; $i<$n; $i++) {
    my $m = $ARGV[$i];
    my $s = 0;
    my $m_a = 0;
    for (my $j=0; $j<$m; $j++) {
      if (defined($r->[$p]) && ($r->[$p] ne "NaN") && ($r->[$p] ne "NA")) {
	$m_a ++;
      }
      $p ++;
    }

    if ($m_a >= 3) {
      $cntval ++;
    }

  }

  if ($cntval == $n) {
    print "$l\n";
  }
  $cnt ++;
}

