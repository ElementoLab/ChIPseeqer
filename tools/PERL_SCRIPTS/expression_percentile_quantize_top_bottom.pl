#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print     Sets::jointab($r);
my $v     = $ta->getColumn(1);
my $p25   = Sets::percentile($v, $ARGV[1]);
my $p75   = Sets::percentile($v, $ARGV[2]);

foreach my $r (@$a_ref) {
  if ($r->[1] <= $p25) {
    print "$r->[0]\t0\n";
  } elsif ($r->[1] < $p75) {
    print "$r->[0]\t1\n";
  } else {
    print "$r->[0]\t2\n";
  }
}

print STDERR "p$ARGV[1] = $p25\n";
print STDERR "p$ARGV[2] = $p75\n";


