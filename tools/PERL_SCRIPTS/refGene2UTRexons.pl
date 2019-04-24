#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;

open IN, $ARGV[0];

my $thisi = undef;
if (defined($ARGV[1]) && ($ARGV[1] ne "")) {
  $thisi = $ARGV[1];
}

while (my $l = <IN>) {
  #print "$l";
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $n    = $a[1];
  my $c    = $a[2];
  my $fr   = $a[3];
  my $e_st = $a[9];
  my $e_en = $a[10];
  my $e_fr = $a[15];

  $e_st =~ s/\,$//;
  $e_en =~ s/\,$//;
  $e_fr =~ s/\,$//;

  my @a_e_st = split /\,/, $e_st;
  my @a_e_en = split /\,/, $e_en;
  my @a_e_fr = split /\,/, $e_fr;

  my $cds_st = $a[6];
  my $cds_en = $a[7];

  my @exons = ();
  for (my $i=0; $i<@a_e_st; $i++) {
    my @a_tmp = ($a_e_st[$i], $a_e_en[$i]);
    push @exons, \@a_tmp;

  }

  # 
  for (my $i=0; $i<@exons; $i++) {
    $exons[$i]->[2] = $exons[$i]->[0];
    $exons[$i]->[3] = $exons[$i]->[1];
  }
 
  for (my $i=0; $i<@exons; $i++) {
    if (($cds_st >= $exons[$i]->[0]) && ($cds_st <= $exons[$i]->[1])) {
      $exons[$i]->[2] = $cds_st;
      $exons[$i]->[3] = $exons[$i]->[1];
      last;
    } else {
      $exons[$i]->[2] = -1;
      $exons[$i]->[3] = -1;
    }	
  }

  for (my $i=@exons-1; $i>=0; $i--) {
    if (($cds_en >= $exons[$i]->[0]) && ($cds_en <= $exons[$i]->[1])) {
      $exons[$i]->[2] = $exons[$i]->[0];
      $exons[$i]->[3] = $cds_en;
      last;
    } else {
      $exons[$i]->[2] = -1;
      $exons[$i]->[3] = -1;
    }	
  }
  
  if ($a[3] eq "-") {
    #@exons  = reverse @exons;
    #@a_e_fr = reverse @a_e_fr;
    if ($ARGV[2] ne "") {
      my $a_ref = Sets::shuffle_array(\@exons);
      @exons = @$a_ref;
    }
    #print STDERR "reevrse\n";
  } elsif ($a[3] ne "+") {
    die "pb\n";
  }

  if ($a[3] eq "+") {
    $a[3] = 1;
  } else {
    $a[3] = -1;
  }
  
  my $i         = undef;  # num exon
  if ($a[3] == 1) {
    $i = 1;
  } else {
    $i = @exons;
  }
  my $j         = 0;
  my $t         = undef;
  my $u3exon_st = undef;
  my $u3exon_en = undef;
  
  foreach my $in (@exons) {
    
    if (($in->[0] != $in->[2]) || ($in->[1] != $in->[3])) {

      if ($in->[0] < $cds_st) {
	$u3exon_st = $in->[0];
	$u3exon_en = Sets::min($in->[1], $cds_st);
	if ($a[3] == 1) {
	  $t         = 5;
	} else {
	  $t         = 3;
	}
	print "$n-$t-$i\t$c\t$u3exon_st\t$u3exon_en\t$a[3]\n";
      } 

      if ($in->[1] > $cds_en) {
	$u3exon_st = Sets::max($cds_en, $in->[0]);
	$u3exon_en = $in->[1];	  
	if ($a[3] == 1) {
	  $t         = 3;
	} else {
	  $t         = 5;
	}
	print "$n-$t-$i\t$c\t$u3exon_st\t$u3exon_en\t$a[3]\n";
      } 

      #else {
      #	die "Problem with + exon\n";
      #}
      
     
      
    }

    if ($a[3] == 1) {
      $i++;
    } else {
      $i--;
    }

    $j ++;

  }
  
}

close IN;

