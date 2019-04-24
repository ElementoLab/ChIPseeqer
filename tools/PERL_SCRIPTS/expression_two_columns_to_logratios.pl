#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
my $n = shift @$r;
print "$n\t" . join("/", @$r); print "\n";

foreach my $r (@$a_ref) {

  next if ( $r->[1]/$r->[2] <= 0);
  my $l = undef;
  if ( abs($r->[2] - $r->[1]) < 0.01) {
    $l = 0; 
  } else {
    if ($ARGV[1] eq "") {
      $l = Sets::log2($r->[1]/$r->[2]);
    } else {
      $l = Sets::log2($r->[2]/$r->[1]);
    }
  }


  my $s = sprintf("%4.3f", $l);
  print "$r->[0]\t$s\n";
}

