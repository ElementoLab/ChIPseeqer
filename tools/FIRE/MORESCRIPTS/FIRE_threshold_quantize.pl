#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {

  print "$r->[0]\t";
  if ($r->[1] > $ARGV[1]) {
    print "1\n";
  } else {
    print "0\n";
  }

}

