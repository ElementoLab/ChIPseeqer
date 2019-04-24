use lib qw(/home/olly/PERL_MODULES);

use Table;

# load mapping
my $t = Table->new;
$t->loadFile($ARGV[0]);
my $a_ref = $t->getArray();


# load translation file 1
$t->loadFile($ARGV[1]);
my $h_ref1 = $t->getIndexKV(0, 1);


# load translation file 2
$t->loadFile($ARGV[2]);
my $h_ref2 = $t->getIndexKV(0, 1);


foreach my $r (@$a_ref) {
    print $h_ref1->{$r->[0]} . "\t" . $h_ref2->{$r->[1]} . "\n";
}
