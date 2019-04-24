#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;

my $a_ref1 = Sets::readSet($ARGV[0]);
my $a_ref2 = Sets::readSet($ARGV[1]);

my $a_union = Sets::getNotIntersection($a_ref1, $a_ref2, 0 );

Sets::printSet($a_union);
