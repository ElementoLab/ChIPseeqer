#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

open OUT1, ">$ARGV[0].exp";
open OUT2, ">$ARGV[0].seq";

print OUT1 "GENE\tEXP\n";

my $i = 0;
foreach my $r (@$a_ref) {

  print OUT1 "seq$i\t$r->[0]\n";
  print OUT2 ">seq$i\n$r->[1]NNNNNNNNN$r->[3]\n\n";
  
  $i ++;
}

close OUT1;
close OUT2;


