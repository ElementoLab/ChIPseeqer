#!/usr/bin/perl

#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;


my $a_ref = Sets::readSet($ARGV[0]);

if ($ARGV[1]) {
    foreach my $r (@$a_ref) {
	$r = abs($r);
    }
}

print Sets::median($a_ref);
print "\n";

