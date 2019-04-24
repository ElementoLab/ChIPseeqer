#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {

  if ($r->[1] <= $ARGV[1]) {
    print "$r->[0]\t0\n";
  } elsif ($r->[1] >= $ARGV[2]) {
    print "$r->[0]\t2\n";
  } else {
    print "$r->[0]\t1\n";
  }

}

