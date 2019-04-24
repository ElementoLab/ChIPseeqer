#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;

# load profieles
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my %M = ();
foreach my $r (@$a_ref) {
  
  $M{$r->[0]}{$r->[1]} = 1;

}


# load expfile

$ta->loadFile($ARGV[0]);
my $a_ref_exp = $ta->getArray();

my $a_ref_e = Sets::shuffle_array($a_ref_exp);

foreach my $m (keys(%M)) {
  print "\t$m";
}
print "\n";
foreach my $r (@$a_ref_e) {
  
  print "$r->[0]";
  
  foreach my $m (keys(%M)) {
    if (defined($M{$m}{$r->[0]})) {
      print "\t1";
    } else {
      print "\t0";
    }
  }

  print "\n";

}
