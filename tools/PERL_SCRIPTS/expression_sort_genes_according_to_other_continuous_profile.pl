#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;

# load profile
$ta->loadFile($ARGV[1]);
$ta->shift;
my $a_ref = $ta->getArray();

my @a_new = sort { $a->[1] <=> $b->[1] } @$a_ref;


# load complete table
$ta->loadFile($ARGV[0]);
my $r = $ta->shift;
print Sets::jointab($r);
my $h_ref_t = $ta->getIndex(0);



foreach my $r (@a_new) {
  if (!defined($h_ref_t->{$r->[0]})) {
    die "$r->[0] not defined.\n";
  }
  print join("\t", @{$h_ref_t->{$r->[0]}}) . "\n";
}


