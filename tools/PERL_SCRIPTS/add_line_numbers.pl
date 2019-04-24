#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";

use Table;
use Sets;

my $cnt = 0;
while ($l1 = <STDIN>) {
    print "$cnt\t$l1";
    $cnt++;
}

