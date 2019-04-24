#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0,1);

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray;


foreach my $r (@$a_ref) {
  print join("\t", @$r) . "\t" . $h_ref->{$r->[0]} . "\n"; 
}

