#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
my $n = shift @$r;
my @a = ($n, 'NAME', 'GWEIGHT', @$r);
print join("\t", @a); print "\n";
print "EWEIGHT\t\t"; print "\t1" x scalar(@$r); print "\n";
foreach my $r (@$a_ref) {
  my $n = shift @$r;
  my @a = ($n, $n, 1, @$r);
  print join("\t", @a); print "\n";
}

