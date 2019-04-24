#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

my $f = shift @ARGV;

my @a_tmp;

open IN, $f;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (Sets::in_array($a[1], @ARGV)) {
    $a[1] = join("/", @ARGV);
  } 
  
  #print join("\t", @a) . "\n";
  push @a_tmp, \@a;

}
close IN;

# for loop
my %H = ();
foreach my $r (@a_tmp) {

  if ((!defined($H{$r->[0]}{$r->[1]}{$r->[2]})) || ($r->[4] > $H{$r->[0]}{$r->[1]}{$r->[2]}->[4]))  {
    $H{$r->[0]}{$r->[1]}{$r->[2]} = $r;
  } 


}

foreach my $k1 (keys(%H)) {

  foreach my $k2 (keys(%{$H{$k1}})) {
    
    foreach my $k3 (keys(%{$H{$k1}{$k2}})) {
      
      print join("\t", @{$H{$k1}{$k2}{$k3}}) . "\n";

    }

  }
  
}
