#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

use refGene;

my $re = refGene->new;
$re->setExt(5000);


my $chr  = "chr3";
my $i1   = 188911274;
my $i2   = 188956636;

$chr = $ARGV[1] if ($ARGV[1] ne "");
$i1  = $ARGV[2] if ($ARGV[2] ne "");
$i2  = $ARGV[3] if ($ARGV[3] ne "");

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  # first locus overlaps with BCL6
  #print "L1 = chr=$a[1]:pos=$a[2]\n";
  #print "L2 = chr=$a[5]:pos=$a[6]\n";
  
  if ($a[1] == 23) {
    $a[1] = "X";
  }
  if ($a[1] == 24) {
    $a[1] = "Y";
  }
  
  if ($a[5] == 23) {
    $a[5] = "X";
  }
  if ($a[5] == 24) {
    $a[5] = "Y";
  }
  
  if ($a[1] ne "0") {
    $a[1] = "chr$a[1]";
  }
  if ($a[5] ne "0") {
    $a[5] = "chr$a[5]";
  }

  if (($a[1] eq $chr) && (Sets::sequencesOverlap($i1, $i2, $a[2], $a[2]+1) > 0)) {

    # gene overlap
    my $a_ref = $re->FindGenesThatOverlapWith($a[5], $a[6], $a[6]+1);
    my @a_txt = ();
    foreach my $r (@$a_ref) {
      push @a_txt, $r->[1] if (!Sets::in_array($r->[1], @a_txt));
    }
    my $txt = join("/", @a_txt);

    print "$a[1]\t$a[2]\t=>\t$a[5]\t$a[6]\t$txt\n";

  } elsif (($a[5] eq $chr) && (Sets::sequencesOverlap($i1, $i2, $a[6], $a[6]+1) > 0)) {
    
    # gene overlap
    my $a_ref = $re->FindGenesThatOverlapWith($a[1], $a[2], $a[2]+1);
    my @a_txt = ();
    foreach my $r (@$a_ref) {
      push @a_txt, $r->[1] if (!Sets::in_array($r->[1], @a_txt));
    }
    my $txt = join("/", @a_txt);

  
    print "$a[5]\t$a[6]\t=>\t$a[1]\t$a[2]\t$txt\n";
  }

}
close IN;

