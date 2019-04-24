#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
my $n = shift @$r;
print "$n\t" . join("_", @$r) . "\n";

foreach my $r (@$a_ref) {

  my $n = shift @$r;

  print "$n";

  my $cntup = 0;
  my $cntdn = 0;
  
  foreach my $s (@$r) {
    if ($s >= $ARGV[1]) {
      $cntup ++;
    } elsif ($s <= 1/$ARGV[1]) {
      $cntdn ++;
    } 
  }

  if ($cntup > @$r/2) {
    print "\t2\n";
  } elsif ($cntdn > @$r/2) {
    print "\t1\n";
  } else {
    print "\t0\n";
  }

}

