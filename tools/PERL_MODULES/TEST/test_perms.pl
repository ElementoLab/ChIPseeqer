BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;

my $a_ref = Sets::get_all_permutations($ARGV[0]);

foreach my $r (@$a_ref) {
  print join(" ", @$r); print "\n";
}
