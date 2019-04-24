#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

system("cp $ARGV[0] $ARGV[1].befrep");
system('perl -pi -e "s/MAL6P1\.287/PFF0670w/g" ' . $ARGV[0]);
system('perl -pi -e "s/MAL6P1\.44/PFF0200c/g" ' . $ARGV[0]);
