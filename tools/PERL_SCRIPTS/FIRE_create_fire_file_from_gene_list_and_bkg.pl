#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;



# load the clusters
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndex(0);


#
# load the background
#
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
print "i\ti\n";
foreach my $r (@$a_ref) {
  if (defined($h_ref->{$r->[0]})) {
    print "$r->[0]\t1\n";
  } else {
    print "$r->[0]\t0\n";
  }
}

