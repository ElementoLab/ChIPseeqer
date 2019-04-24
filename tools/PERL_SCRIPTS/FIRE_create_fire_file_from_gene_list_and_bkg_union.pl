#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;


#
# load the clusters
#
#$ta->loadFile($ARGV[0]);
my $h_ref = Sets::getIndex($ARGV[0]);




#
# load the background
#
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
print "i\ti\n";
foreach my $r (@$a_ref) {
  if (defined($h_ref->{$r->[0]})) {
    print "$r->[0]\t1\n";
    $h_ref->{$r->[0]} = 2; # mark
  } else {
    print "$r->[0]\t0\n";
  }
}

foreach my $k (keys(%$h_ref)) {
  #print "$h_ref->{$k}\n";
  if ($h_ref->{$k} == 1) {
    print "$k\t1\n";
  }
}
