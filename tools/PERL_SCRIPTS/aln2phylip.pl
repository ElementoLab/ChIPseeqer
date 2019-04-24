#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use  ClustalW;
use strict;

my $cl = ClustalW->new;


$cl->setFile($ARGV[0]);

my $hcode = undef;
if ($ARGV[1] ne "") {
  $hcode = {};
}
  

print $cl->getPhylipFormat($hcode);

if (defined($hcode)) {
  open OUT, ">$ARGV[0].phydict" or die "Cannon open $ARGV[0].phydict\n";
  foreach my $k (keys(%$hcode)) {
    print OUT "$k\t$hcode->{$k}\n";
  }
  close OUT;
  print STDERR "$ARGV[0].phydict created !\n";
}
