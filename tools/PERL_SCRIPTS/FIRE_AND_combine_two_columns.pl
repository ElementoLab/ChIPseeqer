#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
#my $r     = shift @$a_ref;
#print Sets::tabjoin($r);

my $i = $ARGV[1];
my $j = $ARGV[2];

foreach my $r (@$a_ref) {

  if (($r->[$i] >= 1) && ($r->[$j] >= 1)) {
    print "$r->[0]\t1\n";
  } else {
    print "$r->[0]\t0\n";
  }

}

