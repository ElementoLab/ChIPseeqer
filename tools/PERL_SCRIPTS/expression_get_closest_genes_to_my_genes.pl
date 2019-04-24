#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndexShifted();


my %H = ();
foreach my $r (keys(%$h_ref)) {

  my $p = 0;
  my $n = 0;
  for (my $i=1; $i<@ARGV; $i++) {
    if (defined($h_ref->{$ARGV[$i]})) {
      $p += Sets::pearson($h_ref->{$ARGV[$i]}, $h_ref->{$r});
      $n++;
    }
  } 

  $H{ $r } = $p/$n;
  print "$r\t$H{$r}\n";
}

#my $o_ref = Sets::hash_order(\%H);
#foreach my $k (reverse(@$o_ref)) {
#  print "$k\t" . sprintf("%4.3f", $H{$k}) . "\n";
#}

