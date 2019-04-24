#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  foreach my $s (@$r) {
    if ($ARGV[1] eq "") {
      $s = "'$s'";
    } else {
      $s = "\"$s\"";
    }
  }
  print join("\t", @$r) . "\n";
}



