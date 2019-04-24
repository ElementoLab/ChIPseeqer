#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";


open IN, $ARGV[0];
my $l = <IN>; chomp $l;
my @a = split /\t/, $l;
foreach my $l (@a) {
  $l =~ s/$ARGV[1].+$//;
}
print join("\t", @a) . "\n";
while (my $l = <IN>) {
  print "$l";
}
close IN;

