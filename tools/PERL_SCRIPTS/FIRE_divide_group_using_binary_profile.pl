#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Hypergeom;
use strict;

my $ta = Table->new;

my $s1 = 0;
my $s2 = 0;
my $ov = 0;
my $nn = 0;

#
# build an index based on the second profile (genes shouls be unique)
#
$ta->loadFile($ARGV[1]);
my $h_ref1 = {};
my $a_ref1 = $ta->getArray();
shift @$a_ref1;
foreach my $r (@$a_ref1) {
  if (defined($h_ref1->{ $r->[0] })) {
    die "Genes in second profile should be unique\n";
  }
  $h_ref1->{ $r->[0] } = $r->[1];
}


$ta->loadFile($ARGV[0]);
my $a_ref2 = $ta->getArray();

my $r = shift @$a_ref2;
print Sets::jointab($r);
foreach my $r (@$a_ref2) {
  if (($r->[1] >= 1) && defined($h_ref1->{$r->[0]})) {
    if ($h_ref1->{$r->[0]} >= 1) {
      print "$r->[0]\t2\n";  # t1 and t2
    } else {
      print "$r->[0]\t1\n";  # t1 and not t2
    }
  } else {
    print "$r->[0]\t0\n";    # not t1
  }
}


