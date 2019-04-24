#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if ($r->[0] ne "VEGA") {

    if (Sets::in_array($r->[6], ("EXP",
				 "IC",
				 "IDA",
				 "IEP",
				 "IGI",
				 "IMP",
				 "IPI",
				 "RCA",
				 "TAS"))) {

      print "$r->[2]\t$r->[4]\n";
    }
  }
}

