BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;


my $a_ref = Sets::readSet($ARGV[0]);


print Sets::stddev($a_ref);
print "\n";

