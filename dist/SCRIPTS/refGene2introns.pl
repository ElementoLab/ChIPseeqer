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

  my $n = $a[1];
  my $c = $a[2];

  my $e_st = $a[9];
  my $e_en = $a[10];
  
  $e_st =~ s/\,$//;
  $e_en =~ s/\,$//;
  
  my @a_e_st = split /\,/, $e_st;
  my @a_e_en = split /\,/, $e_en;

  #my $intron  = 1;
  my @introns = ();
  for (my $i=0; $i<@a_e_en-1; $i++) {
    my @a_tmp = ($a_e_en[$i], $a_e_st[$i+1]);
    push @introns, \@a_tmp;
    #print "$n-$intron\t" . $c . "\t" .  . "\n";
    #$intron ++;
  }

  if ($a[3] eq "-") {
    @introns = reverse @introns;
    if ($ARGV[2] ne "") {
      my $a_ref = Sets::shuffle_array(\@introns);
      @introns = @$a_ref;
    }
    #print STDERR "reevrse\n";
  } elsif ($a[3] ne "+") {
    die "pb\n";
  }

  my $i = 1;
  foreach my $in (@introns) {
    if ((defined($thisi) && (($thisi == 0) || ($thisi == $i))) || !defined($thisi)) { 
      print "$n-$i\t$c\t$in->[0]\t$in->[1]\n";
    }
    
    $i ++;
  }
  
}

close IN;

