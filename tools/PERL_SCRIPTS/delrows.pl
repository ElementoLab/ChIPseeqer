#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;

my $i = 0;
while ($l = <STDIN>) {
    print $l if (!Sets::in_array($i, @ARGV));
    $i ++;
}




