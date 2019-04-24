#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);

use Table;

my $ta = Table->new;

$ta->loadFile($ARGV[0]);

print $ta->getHtmlTable();
