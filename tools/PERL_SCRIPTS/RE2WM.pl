#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;


print "Motif 1\n";
print Sets::myre2scanace($ARGV[0]);
#print "\n";


