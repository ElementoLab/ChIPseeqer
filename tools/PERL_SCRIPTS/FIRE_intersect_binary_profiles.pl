#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";
use lib "/home/ole2001/usr/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";

use Table;
use Hypergeom;
use strict;

my $ta = Table->new;

my $s1 = 0;
my $s2 = 0;
my $ov = 0;
my $nn = 0;

#
# build an index based on the second profile (genes shouls be unique)
#
$ta->loadFile($ARGV[1]);
#$ta->shift;
my $h_ref1 = {};
my $a_ref1 = $ta->getArray();
shift @$a_ref1;
foreach my $r (@$a_ref1) {
  if (defined($h_ref1->{ $r->[0] })) {
    die "Genes in second profile should be unique\n";
  }
  $h_ref1->{ $r->[0] } = $r->[1];
}


$ta->loadFile($ARGV[0]);
$ta->shift;
my $a_ref2 = $ta->getArray();

my @OV = ();
foreach my $r (@$a_ref2) {

  if (defined($h_ref1->{$r->[0]})) {
    $nn ++;
    
    
    if ($r->[1] >= 1) {
      $s1 ++;
    }

    if ($h_ref1->{$r->[0]} >= 1) {
      $s2 ++;
    }

    if (($r->[1] >= 1) && ($h_ref1->{$r->[0]} >= 1)) {
      $ov ++;
      push @OV, $r->[0];
    }
    

  }

}


print "($ov, $s1, $s2, $nn)\n";
print Hypergeom::cumhyper($ov, $s1, $s2, $nn); print "\n";

if ($ARGV[2]) {
  print join("\n", @OV) . "\n";
}
