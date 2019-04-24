#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;
if (@ARGV == 0) {
	die "Args: matrix\n";

}
my $ta = Table->new;

# LOAD PROFILES
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my @SUM = ();

foreach my $r (@$a_ref) {

  my $n = shift @$r;
  for (my $i=0; $i<@$r; $i++) {
    $SUM[$i] += $r->[$i];
  }

  
  
}

print join("\t", @SUM) . "\n";
