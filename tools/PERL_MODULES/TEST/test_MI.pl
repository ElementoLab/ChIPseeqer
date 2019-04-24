#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my @a = (1,0,1,0,0,0,1,1,1,1,0,0,1,0,0);
my @b = (0,0,0,1,0,0,1,1,1,1,0,1,0,1,0);
my @c = (0,0,0,0,0,0,1,1,1,1,0,0,0,0,0);

my @e = (0,0,0,0,0,1,1,1,1,1,2,2,2,2,2);

my $mi = Sets::CalculateMI(\@c, \@e, 2, 3);

print "MI=$mi\n";
