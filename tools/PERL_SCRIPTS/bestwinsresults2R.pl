use lib qw(/home/olly/PERL_MODULES);

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
    my $a= $r->[5];

    my ($x,$y) = $a =~ /\((\d+)\/(\d+)\)/;

    $r->[1] =~ s/\ //g;

    print "$r->[1]\t$x\t$y\n";
	
}
