#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


my $map = "/Users/olivier/PEOPLE/RITA/NORMALS/Olivier/5k_associatedRefSeqUcscGenesHG17toHG19.txt";

my %PR = ();

open IN, $map or die "Cannot open $map\n";
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $ups = $a[11];
  
  next if ($ups =~ /No/);
  next if ($a[4] eq "");

  #push @{$PR{$a[0]}}, $a[3] if !Sets::in_array($a[3], @{$PR{$a[0]}});

  push @{$PR{$a[0]}}, $a[4] if !Sets::in_array($a[3], @{$PR{$a[0]}});
  
}
close IN;

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;
print $l;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1; 
        my $pro = $a[0];
  
  if (defined($PR{$a[0]})) {
    foreach my $g (@{$PR{$a[0]}}) {
      $a[0] = $g;
      print join("\t", @a) . "\t$pro" . "\n";
    }
  }
  
}
close IN;


