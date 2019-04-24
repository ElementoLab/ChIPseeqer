BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {
  push @{ $H{$r->[0]} }, $r->[1];
}


foreach my $n (keys(%H)) {
  print "$n\t";
  print join(",", @{$H{$n}});
  print "\n";
}

