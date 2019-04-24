#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";
if ($ARGV[1] ne "") {
  print STDERR "# setting all values < 1 to 1\n";
} 
foreach my $r (@$a_ref) {
  my $n = shift @$r;
  if ($ARGV[1] ne "") {
    
    foreach my $s (@$r) {
      if ($s < 1) {
        $s = 1;
      }
    }
  }
  $r = Sets::logArray($r);
  print "$n";
  foreach my $s (@$r) {
    print "\t";
    if ($s eq "") {
      print "";
    } else {
      my $prec = 1;
      if ($ENV{PREC} ne "") {
        $prec = $ENV{PREC};
      }
      print sprintf("%3.$prec" . "f", $s);

      #print sprintf("%4.3f", $s);
    }
  }
  print "\n";
}


