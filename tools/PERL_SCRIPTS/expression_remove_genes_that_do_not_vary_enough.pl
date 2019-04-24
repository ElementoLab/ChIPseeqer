#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
if (@ARGV == 0) {
  die "Args: file stdmin\n";
}
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {
  my $m = shift @$r;

  my $s = Sets::stddev($r);

  if ($s >= $ARGV[1]) {
    print "$m\t" . join("\t", @$r) . "\n";
  }
}

