#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

my $np = undef;
open IN, "$ENV{SNVSEEQERDIR}/REFDATA/refLinkrefGene.txt.26Dec2009" or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[0] eq $ARGV[1]) {
    $np = $a[1];
  }
}
close IN;

if (!defined($np)) {
  die "Cannot find $ARGV[1]\n";
}


my @a = `fastacmd -d $ARGV[0] -s $np`;
shift @a;
chomp @a;

my $txt = join("", @a);

my @b = split //, $txt;

print $b[$ARGV[2]-1] . "\n";

