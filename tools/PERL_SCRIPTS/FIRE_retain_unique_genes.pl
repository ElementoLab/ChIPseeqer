#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

my %H = ();
foreach my $r (@$a_ref) {
  if (!defined($H{$r->[0]})) {
    print join("\t", @$r) . "\n";
  } 
  $H{ $r->[0] } = 1;
}

