#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref1 = $ta->getColumn(0);
$ta->loadFile($ARGV[1]);
my $a_ref2 = $ta->getColumn(0);

my $a_union = Sets::getOverlapSet($a_ref1, $a_ref2);

Sets::printSet($a_union);
