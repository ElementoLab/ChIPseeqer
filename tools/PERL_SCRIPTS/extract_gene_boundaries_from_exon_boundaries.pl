use lib qw(/home/elemento/PERL_MODULES);

use Sets;
use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my %H;
foreach my $r (@$a_ref) {
    
    #if (defined($H{$r})) {
    next if (($r->[2] eq "") || ($r->[3] eq ""));
    
    if ($H{$r->[0]}) {
	$H{$r->[0]}->[2] = Sets::min($H{$r->[0]}->[2], $r->[2]);
	$H{$r->[0]}->[3] = Sets::max($H{$r->[0]}->[3], $r->[3]);
    } else {
	my @a = ($r->[0], $r->[1], $r->[2], $r->[3], $r->[4]);
	$H{$r->[0]} = \@a;
    }

}

foreach my $r (sort(keys(%H))) {
    print join("\t", @{ $H{$r} });
    print "\n";
}


