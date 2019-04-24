#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print     Sets::jointab($r);
#print "GENE\tEXP\n";

foreach my $r (@$a_ref) {
  if ($r->[1] >= 1) {
    $r->[1] = 1
  }
  print "$r->[0]\t$r->[1]\n";
}

