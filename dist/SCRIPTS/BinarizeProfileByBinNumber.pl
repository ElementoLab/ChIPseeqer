#!/usr/bin/perl
use lib "$ENV{CHIPSEEQERDIR}";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;

print join("\t", @$r);
print "\n";

my @a_bins = split /\,/, $ARGV[1];

foreach my $r (@$a_ref) {
  
  if (Sets::in_array($r->[1], @a_bins)) {
    $r->[1] = 1;
  } else {
    $r->[1] = 0;
  }

  print join("\t", @$r);
  print "\n";
}

