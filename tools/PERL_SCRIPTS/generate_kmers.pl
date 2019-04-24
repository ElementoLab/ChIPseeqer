#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;


#Sets::printSet(Sets::allkmers($ARGV[0]));

#Sets::printSet(

#my $a_ref = Sets::removeLowComplexityKmers(, $ARGV[0]-1);


my $a_ref1 = Sets::allkmers($ARGV[0]);
if (defined($ARGV[1])) {
  $a_ref1 = Sets::removeComplements($a_ref1);
}

Sets::printSet($a_ref1);



