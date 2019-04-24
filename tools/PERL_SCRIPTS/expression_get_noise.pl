#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;

foreach my $r (@$a_ref) {
  my $n  = shift @$r;
  my $ra = Sets::stddev($r) / Sets::average($r);
  print "$n\t$ra\n";
}

