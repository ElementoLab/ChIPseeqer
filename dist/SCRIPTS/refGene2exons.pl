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
    @exons  = reverse @exons;
    @a_e_fr = reverse @a_e_fr;
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
  
  my $i = 1;
  my $j = 0;
  foreach my $in (@exons) {
    if ((defined($thisi) && (($thisi == 0) || ($thisi == $i))) || !defined($thisi)) { 
      print "$n-$i\t$c\t$in->[0]\t$in->[1]\t$a[3]\t$a_e_fr[$j]\t$in->[2]\t$in->[3]\n";
    }
    $j ++;
    $i ++;
  }
  
}

close IN;

