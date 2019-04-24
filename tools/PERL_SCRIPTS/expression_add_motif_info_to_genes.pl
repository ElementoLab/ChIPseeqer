#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();


my %HASM = ();
foreach my $r (@$a_ref) {
  $HASM{ $r->[0] } .= $r->[1];
}



$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  $r->[0] = "$HASM{$r->[0]} $r->[0]";
  print join("\t", @$r) . "\n";
}


