#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my @sums = ();
my @nts  = ();
foreach my $r (@$a_ref) {
  my $n = shift @$r;
  push @nts, $n;
  for (my $i=0; $i<@$r; $i++) {
    $sums[$i] += $r->[$i];
  }
}


#print join("\t", @sums);
#print "\n";

my $j = 0;
foreach my $r (@$a_ref) {
  #my $n = shift @$r;
  #push @nts, $n;
  print "$nts[$j]";
  for (my $i=0; $i<@$r; $i++) {
    print "\t". sprintf("%4.3f", $r->[$i]/$sums[$i]);
  }
  print "\n";
  $j++
}


