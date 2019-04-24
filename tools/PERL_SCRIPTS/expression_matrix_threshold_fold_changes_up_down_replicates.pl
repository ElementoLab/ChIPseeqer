#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
my $n = shift @$r;
print "$n\t" . join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {

  my $n = shift @$r;

  my $cntup = 0;
  my $cntdn = 0;
  
  foreach my $s (@$r) {
    if (($ARGV[1] > 0) && ($s >= $ARGV[1])) {
      $cntup ++;
    } elsif (($ARGV[1] < 0) && ($s <= -1/$ARGV[1])) {
      $cntdn ++;
    } 
  }

  if ($cntup > @$r/2) {
    print "$n\t" . join("\t", @$r) . "\n";
  } elsif ($cntdn > @$r/2) {
    print "$n\t" . join("\t", @$r) . "\n";
  } else {
    #print "\t0\n";
  }

}

