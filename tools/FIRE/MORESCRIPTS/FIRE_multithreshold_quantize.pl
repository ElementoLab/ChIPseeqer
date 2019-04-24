#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

shift @ARGV;
my @t = @ARGV;

foreach my $r (@$a_ref) {

  print "$r->[0]\t";

  if ($r->[1] < $t[0]) {
    print "0\n";
  } else {
    for (my $i=0; $i<@t; $i++) {
      if ((($i == scalar(@t)-1) && ($r->[1] >= $t[$i])) ||
	  (($r->[1] >= $t[$i]) && ($r->[1] < $t[$i+1]))) {
	my $ii = $i+1;
	
	print "$ii\n";
	last;
      } 
    }
  }


}

