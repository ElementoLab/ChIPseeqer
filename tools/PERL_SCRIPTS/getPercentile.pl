#!/usr/bin/perl

#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;


my $a_ref = Sets::readSet($ARGV[0]);

print Sets::percentile($a_ref, $ARGV[1]);
print "\n";

