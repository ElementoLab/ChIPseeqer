#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;

my $f = shift @ARGV;

$ta->loadFile($f);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";
foreach my $r (@$a_ref) {
  print join("\t", @$r) . "\n" if (Sets::in_array($r->[1], @ARGV));
}

