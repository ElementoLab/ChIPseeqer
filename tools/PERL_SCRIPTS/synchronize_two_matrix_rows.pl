#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";
use Table;
use Sets;
use strict;

if (@ARGV == 0) {
  die("Args: mat1 mat2\n");
}

my $ta = Table->new;
$ta->loadFile($ARGV[0]) or "Cannot load matrix1 file\n";
my $a_ref1_h = $ta->shift;
my $a_ref1_c = $ta->getColumn(0);
my $h_ref1   = $ta->getIndex(0);

$ta->loadFile($ARGV[1]) or "Cannot load matrix2 file\n";
my $a_ref2_h = $ta->shift;
my $a_ref2_c = $ta->getColumn(0);
my $h_ref2   = $ta->getIndex(0);


my $ov = Sets::getOverlapSet($a_ref1_c, $a_ref2_c);

my $ov1 = Sets::getNotIntersection($a_ref1_c, $a_ref2_c);
print STDERR join(" ", sort(@$ov1)) . "\n";

my $ov2 = Sets::getNotIntersection($a_ref2_c, $a_ref1_c);
print STDERR join(" ", sort(@$ov2)) . "\n";

my $out1 = Sets::filename($ARGV[0]) . ".sync";
my $out2 = Sets::filename($ARGV[1]) . ".sync";

if (-e $out1) {
  print "$out1 already exists, overwrite?\n";
  my $ans = <STDIN>;
  if ($ans =~ /n/) {
    exit;
  }
}

if (-e $out2) {
  print "$out2 already exists, overwrite?\n";
  my $ans = <STDIN>;
  if ($ans =~ /n/) {
    exit;
  }
}

open OUT1, ">$out1";
print OUT1 join("\t", @$a_ref1_h) . "\n";
open OUT2, ">$out2";
print OUT2 join("\t", @$a_ref2_h) . "\n";

foreach my $r (@$ov) {
  print OUT1 join("\t", @{$h_ref1->{$r}}) . "\n";
  print OUT2 join("\t", @{$h_ref2->{$r}}) . "\n";
}

close OUT1;
close OUT2;


if (-e $out1) {
  print "Created $out1\n";
}
if (-e $out2) {
  print "Created $out2\n";
}
