#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getColumn(0);

$ta->loadFile($ARGV[1]);
my $a_ref_m = $ta->getArray();
my $r = shift @$a_ref_m;
my $h_ref   = $ta->getIndex(0);

print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {
  
  if (defined($h_ref->{$r})) {
    print join("\t", @{$h_ref->{$r}}) . "\n";
  }

}

