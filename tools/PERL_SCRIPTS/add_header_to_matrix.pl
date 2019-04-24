use lib qw(/home/elemento/PERL_MODULES);
use Sets;

my $a_ref = Sets::readSet($ARGV[0]);

print "\t" . join("\t", @$a_ref); print "\n";
