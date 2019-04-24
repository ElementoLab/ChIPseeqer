use lib qw(/home/olly/PERL_MODULES);

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my %H = ();
foreach my $r (@$a_ref) {
    
    my @b = split /\//, $r->[1], -1;

    foreach my $s (@b) {
	
	push @{ $H{ $s } }, $r->[0] if (!Sets::in_array($r->[0], @{ $H{ $s } }));

    }

    #print join("\/", );
}


foreach my $k (keys(%H)) {
    print "$k\t"; print join("\t", @{ $H{ $k } }); print "\n";
}
