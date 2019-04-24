#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $a_ref = Sets::readSet($ARGV[1]);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndexKV(0,1);

print "GENE\tEXP\n";
foreach my $r (@$a_ref) {
  if (defined($h_ref->{$r})) {
    print "$r\t" . $h_ref->{$r} . "\n";
  } else {
    print "$r\t0\n";
  }
}

