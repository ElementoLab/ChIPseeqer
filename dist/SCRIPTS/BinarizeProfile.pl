#!/usr/bin/perl
use lib "$ENV{CHIPSEEQERDIR}";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print     Sets::jointab($r);
#print "GENE\tEXP\n";

foreach my $r (@$a_ref) {
  if ($r->[1] >= 1) {
    $r->[1] = 1;
  }
  print "$r->[0]\t$r->[1]\n";
}

