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
my %WM  = ();
foreach my $r (@$a_ref) {

  my $n = shift @$r;
  for (my $i=0; $i<@$r; $i++) {
    $SUM[$i] += $r->[$i];
  }
  $WM{$n} = $r;

}

print "Protein: $ARGV[0]\n";

foreach my $a (("A", "C", "G", "T")) {
  for (my $i=0; $i<@{$WM{$a}}; $i++) {
    $WM{$a}->[$i] /= $SUM[$i];
    $WM{$a}->[$i]  = sprintf("%4.3f", $WM{$a}->[$i]);
  }
  print "$a:\t" . join("\t", @{$WM{$a}}) . "\n";
}

#print join("\t", @SUM) . "\n";
