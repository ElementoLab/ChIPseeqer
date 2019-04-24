#!/usr/bin/perl
use lib "$ENV{CHIPSEEQERDIR}";

use Table;

my $ta = Table->new;

# load first file
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

# load second file
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0,1);

shift @$a_ref;
print "GENE\tEXP\n";
foreach my $r (@$a_ref) {

  next if (!defined($h_ref->{ $r->[0] }));

  my $i = 0;
  if (($r->[1] == 0) && ($h_ref->{ $r->[0] } == 0)) {
    $i = 0;
  } elsif (($r->[1] == 0) && ($h_ref->{ $r->[0] } == 1)) {
    $i = 0;
  } elsif (($r->[1] == 1) && ($h_ref->{ $r->[0] } == 0)) {
    $i = 1;
  } elsif (($r->[1] == 1) && ($h_ref->{ $r->[0] } == 1)) {
    $i = 0;
  }
  print "$r->[0]\t$i\n";

}
