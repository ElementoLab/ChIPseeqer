#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use strict;

use Sets;

my %ROWS = ();


open IN, $ARGV[0];
my $l = <IN>;
print $l;
my $cnt = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{ $ROWS{ $a[0] } }, $a[1];
  $cnt ++;
}

close IN;

#
# 
#
my @LL = sort(keys(%ROWS));

foreach my $k (@LL) {
  my $r = $ROWS{$k};
  my $m = Sets::maxInArray($r);
  print "$k\t$m\n";
}
