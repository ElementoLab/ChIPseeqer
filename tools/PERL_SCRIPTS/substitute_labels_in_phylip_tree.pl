#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";
use Table;
#
#  in: phylip tree, labels, out : phylip tree
#


open IN, "$ARGV[0]";
my @a = <IN>; 
my $l = join("", @a);
close IN;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {

    my $s1 = $r->[0];
    my $s2 = $r->[1];

    $l =~ s/$s1\:/$s2\:/g;
}


print $l;
