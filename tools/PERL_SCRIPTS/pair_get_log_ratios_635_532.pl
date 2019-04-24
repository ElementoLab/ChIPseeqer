#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use strict;
use Sets;

my %V0 = ();

open IN, $ARGV[0];
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{$V0{$a[2]}}, $a[9]; 
}
close IN;

my %V1 = ();
open IN, $ARGV[1];
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{$V1{$a[2]}}, $a[9]; 
  #if (abs($V1{$a[3]}) < 0.001) {
  #  $V1{$a[3]} = 0.001;
  #}
}
close IN;

print "PROBESET\tEXP\n";
foreach my $g (keys(%V0)) {
  next if ($g =~ /RANDOM/);
  my $av = Sets::logratioArrayEntries( $V0{$g}, $V1{$g} );
  my $l  = Sets::average($av);
  #my $m0 = Sets::median($V0{$g});
  #my $m1 = Sets::median($V1{$g});
  #  my $l  = log( $m0/$m1 );
  print "$g\t$l\n"; #($V0{$g}/$V1{$g})\n";


}

