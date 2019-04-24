use lib qw(/home/olly/PERL_MODULES);

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my %H;
foreach my $r (@$a_ref) {
    
    my $len = undef;
    if ($r->[4] == 1) {
	$len = $r->[6] - $r->[3];
    } else {
	$len = $r->[2] - $r->[5];
    }
    
    print "$r->[0]\t$len\n";

}
