#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use lib "$ENV{HOME}/usr/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";

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
$ta->loadFile($ARGV[0]);
$ta->shift;
my $a_ref2 = $ta->getArray();

my @OV = ();
foreach my $r (@$a_ref2) {

  $nn ++;

  if ($r->[$ARGV[1]] >= 1) {
    $s1 ++;
  }

  if ($r->[$ARGV[2]] >= 1) {
    $s2 ++;
  }

  if (($r->[$ARGV[1]] >= 1) && ($r->[$ARGV[2]] >= 1)) {
    $ov ++;
    push @OV, $r->[0];
  }

}


print "($ov, $s1, $s2, $nn)\n";
print Hypergeom::cumhyper($ov, $s1, $s2, $nn); print "\n";

if ($ARGV[3]) {
  print join("\n", @OV) . "\n";
}
