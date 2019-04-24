#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;


my $ta = Table->new;

#
# load matrix from which to get gene
#
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndex(0);

#
# load matrix to display
#
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r = shift @$a_ref;

print join("\t", @$r) . "\n";

my $r = $h_ref->{ $ARGV[2] };
#for (my $i=1; $i<@$r; $i++) {
#$r->[$i] = $r->[$i] / 4;
#}
print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {
  print join("\t", @$r) . "\n";
}

