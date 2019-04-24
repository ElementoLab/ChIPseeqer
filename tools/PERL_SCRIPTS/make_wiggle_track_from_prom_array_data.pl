#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

if (@ARGV == 0) {
  die "Args: suffix coords desc\n";
}

#
# read data
#
#my $ta = Table->new;
#$ta->loadFile($ARGV[0]);
#my $h_ref = $ta->getIndexKV(0,1);



my $suf = $ARGV[0];

#
# load 532 data
#
my $f532 = "$suf\_532.pair";
if (! -e $f532) {
  $f532 .= ".txt";
}

#
# load 635 data
#
my $f635 = "$suf\_635.pair";
if (! -e $f635) {
  $f635 .= ".txt";
}

my %H1 = ();
open IN, $f635;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  if ($a[9] == 0.0) {
    $a[9] = 0.1;
  }

  $H1{$a[3]} = $a[9];
}
close IN;

my $h_ref = {};

open IN, $f532;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;

  #print "$a[3]\t" . $H1{$a[3]} . "\t$a[9]\n";

  if ($a[9] == 0.0) {
    $a[9] = 0.1;
  }

  my $rat = sprintf("%4.3f", log($H1{$a[3]} / $a[9]));
  $h_ref->{$a[3]} = $rat;
  
}
close IN;



# read and store coords
my %CHR = ();
open IN, $ARGV[1];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{ $CHR{ $a[1] } }, \@a;
}
close IN;

# sort per chr
#print "track type=wiggle_0 name=\"$c probes\" description=\"$c probes\"\n";
print "track type=wiggle_0 name=\"$ARGV[2]\" description=\"$ARGV[2]\"\n";
foreach my $c (keys(%CHR)) {
  

  print "variableStep chrom=$c span=50\n";

  my $a_ref = $CHR{$c};
  my @a_sor = sort { $a->[2] <=> $b->[2] } @$a_ref;

  
  my $prev = undef;
  foreach my $p (@a_sor) {
    next if ($p->[2] == $prev);
    my $v = $h_ref->{$p->[0]};
    $v = exp($v);
    print "$p->[2]\t$v\n";    #join("\t", @$p) . "\n";
    $prev = $p->[2];
  }


}





