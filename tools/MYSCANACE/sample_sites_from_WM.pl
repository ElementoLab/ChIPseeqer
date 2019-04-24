#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Motif;

if (scalar(@ARGV) == 0) {
  die "Usage: perl sample_sites_from_WM.pl wmfile\n";
}

my $mo = Motif->new;

my $n  = (defined($ARGV[1])?$ARGV[1]:10000);

$mo->readClassicalWM($ARGV[0]);
# $mo->print;

print "Motif 1\n";
for (my $i=0; $i<$n; $i++) {
  my $s = $mo->generateRandomSequence();
  print "$s\n";
}
for (my $j=0; $j<$mo->getSize(); $j++) {
  print "*";
}
print "\n";
