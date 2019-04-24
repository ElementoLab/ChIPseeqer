#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my $st  = $r->[4];
  my $tss = undef;
  if ($st == 1) {
    $tss = $r->[5];
  } else {
    $tss = $r->[6];
  }
  print "$r->[0]\t$r->[1]\t$st\t$tss\n";
}

