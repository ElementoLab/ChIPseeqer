#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

die "Define a number of genes.\n" if ($ARGV[1] eq "");

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
$ta->shift;
$ta->sortbycol(1);

my $a_ref = $ta->getArray();

print "GENE\tEXP\n";
my $i = 0;
foreach my $r (@$a_ref) {
  if ($i < $ARGV[1]) {
    $r->[1] = 1;
  } else {
    $r->[1] = 0;
  }
  print join("\t", @$r) . "\n";

  $i++;
}

