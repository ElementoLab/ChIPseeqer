#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;

print Sets::jointab($r);

foreach my $r (@$a_ref) {
  $r->[1] = -1 * $r->[1];
  print Sets::jointab($r);
}

