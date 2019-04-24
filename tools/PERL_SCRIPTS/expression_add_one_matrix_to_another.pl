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
shift @$h;

open IN, $ARGV[0];

# first row is special
my $l = <IN>; chomp $l;
my @b = split /\t/, $l, -1;
shift @b;

print "GENE\t" . join("\t", @b) . "\t" . join("\t", @$h) . "\n";
my $m = @$h;

while (my $l = <IN>) {

  chomp $l;
  my @a = split /\t/, $l, -1;
  my $n = shift @a;
  print "$n\t" .  join("\t", @a);

  if (defined($h_ref->{$n})) {
    print "\t" . join("\t", @{$h_ref->{$n}});
  } else {
    for (my $i=0; $i<$m; $i++) {
      print "\t";
    }
  }

  print "\n";
  
}

close IN;

