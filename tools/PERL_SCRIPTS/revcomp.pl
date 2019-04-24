#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";



use Sets;


print Sets::getComplement($ARGV[0]);
print "\n";


