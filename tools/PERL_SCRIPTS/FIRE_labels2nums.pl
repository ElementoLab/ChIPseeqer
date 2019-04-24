#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

# first pass to collect labels
my %LABELS = ();
open IN, $ARGV[0] or die "Cannot open file";
my $l =  <IN>;
#print "$l\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $LABELS{$a[1]} = 1;
}
close IN;

my @lab = sort(keys(%LABELS));
my %H = ();
my $i = 0;
open OUT, ">$ARGV[0].LABELS";
foreach my $k (@lab) {
  $H{$k} = $i;
  print OUT "$i\t$k\n";
  $i++;
}
close OUT;
if ( -e "$ARGV[0].LABELS") {
  print STDERR "Wrote $ARGV[0].LABELS\n";
}


open IN, $ARGV[0] or die "Cannot open file";
my $l =  <IN>;
print "$l\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  print "$a[0]\t$H{$a[1]}\n";

}
close IN;


