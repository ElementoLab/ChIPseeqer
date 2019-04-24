#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;
use Table;
use strict;

# read in the first matrix
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
$ta->processHeader();
my $h     = $ta->getHeader();
my $h_ref = $ta->getIndexShifted();



open IN, $ARGV[0];

# first row is special
my $l = <IN>; chomp $l;
my @b = split /\t/, $l, -1;
shift @b;

print join("\t", @$h) . "\t" . join("\t", @b) . "\n";

while (my $l = <IN>) {

  chomp $l;
  my @a = split /\t/, $l, -1;
  my $n = shift @a;
  if (defined($h_ref->{$n})) {
    print "$n\t" . join("\t", @{$h_ref->{$n}}) . "\t" . join("\t", @a) . "\n";
  }

}

close IN;




