BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my @MOT = ();

foreach my $r (@$a_ref) {
  if ($r->[5] ne "OK-NO-SEED") {
    my @a = ($r->[0], $r->[1], 0);
    push @MOT, \@a;
  }
}


$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if ($r->[5] ne "OK-NO-SEED") {
    my @a = ($r->[0], $r->[1], 1);
    push @MOT, \@a;
  }
}


@MOT = sort { $b->[1] <=> $a->[1] } @MOT;
foreach my $r (@MOT) {
  print "$r->[0]\t$r->[1]\t$r->[2]\n";
}
