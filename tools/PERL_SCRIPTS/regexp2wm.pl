#!/usr/bin/perl
use strict;

if ($ARGV[1] ne "") {
  print myre2ace($ARGV[0]);
} else {
  print myre2wm($ARGV[0]);
}


sub myre2wm {
  my ($re) = @_;

  my @a = split //, $re;
  
  my %WMH = ();
  
  my $i = 0; my $n = length($re);

  my $ml = 0;
  while ($i < $n) {

    # N
    if ($a[$i] eq '.') {
      $WMH{"A"}[$ml] = 1/4;
      $WMH{"C"}[$ml] = 1/4;
      $WMH{"G"}[$ml] = 1/4;
      $WMH{"T"}[$ml] = 1/4;
      $ml++;
    }
    
    # [..]
    elsif ($a[$i] eq '[') {
      
      
      my @b = ();
      $i++;
      while ($a[$i] ne ']') {
        push @b, $a[$i];
        $i ++; 
      }
      foreach my $bb (@b) {
	$WMH{$bb}[$ml] = 1/@b;
      }
      $ml++;

    } else {
      $WMH{$a[$i]}[$ml] = 1.0;
      $ml++;
    }
    
    $i ++;
  }	

  my $txt = "";
  foreach my $nt (("A", "C", "G", "T")) {
    my @ws = ();
    $txt .= "$nt";
    for (my $i=0; $i<$ml; $i++) {
      my $w = $WMH{$nt}[$i];
      if (!defined($w)) {
	$w = 0.0;
      }
      $w = sprintf("%3.2f", $w);
      $txt .= "\t$w";
    }
    $txt .= "\n";
  }
  return $txt;
}



sub myre2ace {
  my ($re) = @_;

  my @a = split //, $re;
  
  my @lines = ();
  
  my $i = 0; my $n = length($re);

  my $ml = 0;
  while ($i < $n) {
    
    
    # N
    if ($a[$i] eq '.') {
      my $l = "";
      $l .= 'A' x 6;
      $l .= 'C' x 6;
      $l .= 'G' x 6;
      $l .= 'T' x 6;
      push @lines, $l;
      $ml++;
    }
    
    # [..]
    elsif ($a[$i] eq '[') {
      
      
      my @b = ();
      $i++;
      while ($a[$i] ne ']') {
        push @b, $a[$i];
        $i ++; 
      }
      
      
      my $m = scalar( @b );
      if ($m == 2) {
        my $l = "";
        foreach my $s (@b) {
          $l .= $s x 12;
        }
        push @lines, $l;
	
      } elsif ($m == 3) {
	my $l = "";
	$l .= $b[0] x 8;
	$l .= $b[1] x 8;
	$l .= $b[2] x 8;
	push @lines, $l;
      }
      $ml++;
    } else {
      my $l = $a[$i] x 24;
      push @lines, $l;
      $ml++;
    }
    
    $i ++;
  }	
  
  #print "ml=$ml\n";
  
  my @M = ();
  my $col = 0;
  foreach my $l (@lines) {
    #print "$l\n";
    my @a = split //, $l; 
    for (my $i=0; $i<24; $i++) {
      $M[ $i ] [ $col ] = $a[$i];
    }
    $col ++;
  } 
  
  my $mo = "";
  for ( my $i=0; $i<24; $i++) {
    for (my $j=0; $j<$ml; $j++) {
      $mo .= $M[$i][$j]; 
    }
    $mo .= "\n";
  }
  
  return $mo;
}


