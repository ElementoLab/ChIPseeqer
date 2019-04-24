#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


# read refgene, gather isoforms
open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";

my %ISOFORMS = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $nm   = $a[1];
  my $ch   = $a[2];
  my $t_st = $a[4] + 1;
  my $t_en = $a[5];
  my $ge   = $a[12];

  my $coords = "$ge\t$ch\t$t_st\t$t_en";
  push @{$ISOFORMS{$coords}}, $nm;

}
close IN;

# process isoforms, ie invert
my %NM2PRI = ();
foreach my $c (keys(%ISOFORMS)) {
  my @a = split /\t/, $c;
  my $idx = 1;
  #print "$c => " . join("\t", @{$ISOFORMS{$c}}) . "\n";
  foreach my $t (@{$ISOFORMS{$c}}) {
    push @{$NM2PRI{$t}}, ["$a[0].$idx", $a[1], $a[2], $a[3] ];  # complete transcript desc
    #print "Adding [$a[0].$idx, $a[1], $a[2], $a[3] ]\n";
    $idx ++;
  }
}




open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";



while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  #next if ($a[2] ne "exon");

  my $e_ch = $a[0];
  my $e_st = $a[3];
  my $e_en = $a[4];
  my $e_fr = $a[6];
  
  my $desc = $a[8];

  my ($gene, $ts) = $desc =~ /gene_id \"(.+?)\"; transcript_id \"(.+?)\"/;

  # determine exact transcript (ts id)
  my $idx = 0;
  my $ovidx = -1;
  foreach my $s (@{$NM2PRI{$ts}}) {
    next if ($s->[1] ne $e_ch);
    if (Sets::sequencesOverlap($e_st, $e_en, $s->[2], $s->[3])) {
      $ovidx = $idx;
    }
    $idx ++;
  }
  my $tss_id = "$ts";
  if ($ovidx == -1) {
    print STDERR  "# Problem with $l\n";
  } else { 
    $tss_id = $NM2PRI{$ts}->[$ovidx]->[0];    
  }
  
  print "$l" . " tss_id \"$tss_id\";\n";

}
close IN;

