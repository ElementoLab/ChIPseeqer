use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getColumn(0);


$ta->loadFile($ARGV[1]);
my $a_ref_gapped   = $ta->getColumn(0);

foreach my $r (@$a_ref_gapped) {
    $r =~ s/N/\./g;
}

my $i = 1;
foreach my $r (@$a_ref_gapped) {

    # nb of matching gapped kmers
    my @a = grep(/$r/, @$a_ref);

    if (scalar(@a) == 0) {
	print "$i\t$r\n";
    }
    $i++;
}
