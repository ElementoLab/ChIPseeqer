use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my $i = 1;
foreach my $r (@$a_ref) {
    $r->[0] = lc($r->[0]);
    
    print ">SEQ$i\n$r->[0]\n";
    $i++;

}
