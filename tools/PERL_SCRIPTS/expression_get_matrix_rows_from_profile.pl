#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

my $prof = shift @ARGV;
my $expf = shift @ARGV;

my $ta = Table->new;
$ta->loadFile($prof);
my $b_ref = $ta->getArray();

my $a_ref = [];
my $h_ref = {};
foreach my $r (@$b_ref) {
  if (Sets::in_array($r->[1], @ARGV)) {
    push @$a_ref, $r->[0];
    $h_ref->{$r->[0]} = $r->[1];
  }
}



#my $a_ref = Sets::readSet($ARGV[0]);
#my $h_ref = Sets::getIndex($ARGV[0]);

my %H = ();

open IN, $expf or die "Cnnot pen $expf\n";
my $m = undef;

my $l = <IN>;
print "$l";
chomp $l;
my @a  = split /\t/, $l, -1;
$m = @a - 1;		


while (my $l = <IN>) {

  chomp $l;
  my @a = split /\t/, $l, -1;

  if (Sets::in_array($h_ref->{ $a[0] }, @ARGV)) { 
    #print "$l\n";
    $H{ $a[0] } = $l;
  } else {
    #print STDERR "$a[0] not found.\n";
  }

}
close IN;

foreach my $r (@$a_ref) {

  if (defined($H{$r})) {
      print "$H{$r}\n";
    
  } else {
    print STDERR "WARNING: $r not found\n";
    #print "$r";
    #  for (my $i=0; $i<$m; $i++) {
    #    print "\t";
    #  }
    #  print "\n";

  }
}
