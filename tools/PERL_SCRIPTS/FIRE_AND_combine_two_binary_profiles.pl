#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0,1);

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;

print Sets::tabjoin($r);

foreach my $r (@$a_ref) {
  if (($r->[1] >= 1) && ($h_ref->{ $r->[0] } >= 1)) {
    print "$r->[0]\t1\n";
  } else {
    print "$r->[0]\t0\n";
  }

}

