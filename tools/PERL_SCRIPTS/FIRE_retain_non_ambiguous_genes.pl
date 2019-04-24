#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;


print join("\t", @$r) . "\n";

my %H = ();
foreach my $r (@$a_ref) {
  $H{$r->[0]} ++;
}

foreach my $r (@$a_ref) {
  if ($H{$r->[0]} == 1) {
    print join("\t", @$r) . "\n";
  }
}

