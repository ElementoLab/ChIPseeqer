#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $b_ref = $ta->getArray();

my $a_ref = [];
my $h_ref = {};
my $h_ref_i = {};
foreach my $r (@$b_ref) {
  push @$a_ref, $r->[0];
  $h_ref->{$r->[0]} = 1;
  $h_ref_i->{$r->[0]} = $r->[1];
}



#my $a_ref = Sets::readSet($ARGV[0]);
#my $h_ref = Sets::getIndex($ARGV[0]);

my %H = ();

open IN, $ARGV[1];
my $m = undef;
if (!defined($ARGV[2])) {
  my $l = <IN>;
  print "$l";
  chomp $l;
  my @a  = split /\t/, $l, -1;
  $m = @a - 1;		
}

while (my $l = <IN>) {

  chomp $l;
  my @a = split /\t/, $l, -1;

  if ($h_ref->{ $a[0] } == 1) {
    #print "$l\n";
    $H{ $a[0] } = $l;
  } else {
    #print STDERR "$a[0] not found.\n";
  }

}
close IN;

foreach my $r (@$a_ref) {

  if (defined($H{$r})) {
    if (defined($h_ref_i->{$r}) && ($h_ref_i->{$r} ne "")) {
      my @a = split /\t/, $H{$r};
      #$a[0] .= "/$h_ref_i->{$r}";
      print join("\t", @a) . "\n";
    } else {
      print "$H{$r}\n";
    }
  } else {
    print STDERR "WARNING: $r not found\n";
    print "$r";
    for (my $i=0; $i<$m; $i++) {
      print "\t";
    }
    print "\n";

  }
}
