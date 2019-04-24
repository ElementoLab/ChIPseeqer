#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my $col = 1;
if (defined($ARGV[1])) {
  $col = $ARGV[1];
}

foreach my $r (@$a_ref) {
  if ($r->[$col] ne "") {
    print join("\t", @$r) . "\n";
  }
}

