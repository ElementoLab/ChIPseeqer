#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";


use Sets;

my $txt = Sets::consensus2re($ARGV[0]);

print "$txt\n";
