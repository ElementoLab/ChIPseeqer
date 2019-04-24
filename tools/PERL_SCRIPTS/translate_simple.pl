use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndexKV(0,1);

while (my $l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;

    if (defined($h_ref->{ $a[0] })) {
	$a[0] = $h_ref->{ $a[0] }
    } else {
	$a[0] = "NA";
    }
    print join("\t", @a); print "\n";
}




