use lib qw(/home/elemento/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref1 = $ta->getArray();

$ta->loadFile($ARGV[1]);
my $h_ref2 = $ta->getIndex(0);

my $max = 100;

foreach my $r (@$a_ref1) {
  if (defined($h_ref2->{ $r->[0] })) {
    print "$r->[0] should be $r->[2], $r->[3], is $h_ref2->{$r->[0]}->[2], $h_ref2->{$r->[0]}->[3]\n";
    if ( (abs($r->[2] - $h_ref2->{$r->[0]}->[2]) <= $max ) && ( abs($r->[3] - $h_ref2->{$r->[0]}->[3] ) <= $max )) {
      print "OK\n";
    }
  }
}
