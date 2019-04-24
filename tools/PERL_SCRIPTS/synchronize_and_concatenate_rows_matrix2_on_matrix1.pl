#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);
use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]) or "Cannot load matrix1 file\n";
my $h_ref1 = $ta->getIndex(0);
my $a_ref1 = $ta->getColumn(0);

$ta->loadFile($ARGV[1]) or "Cannot load matrix2 file\n";
my $h_ref2 = $ta->getIndex(0);

foreach my $r (@$a_ref1) {
  if (defined($h_ref2->{ $r })) {
    print join("\t", @{ $h_ref1->{ $r } }); print "\t";
    shift @{ $h_ref2->{ $r } };
    print join("\t", @{ $h_ref2->{ $r } }); print "\n";
  } else {

    print STDERR "Cannot find $r in matrix 2\n";
  }
}


