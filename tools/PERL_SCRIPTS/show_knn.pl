use lib qw(/home/elemento/PERL_MODULES);
use Sets;
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $h_ref = $ta->getIndex(0);

my $a_ref_row1 = $h_ref->{ $ARGV[1] }; 
shift @$a_ref_row1;


my @SORT = ();
foreach my $r (@$a_ref) {

  my $n = shift @$r;

  my $d = Sets::pearson($a_ref_row1, $r);

  my @a_tmp = ($n, $d);

  push @SORT, \@a_tmp;
  
}


@SORT = sort { $a->[1] <=> $b->[1] } @SORT;

foreach my $s (@SORT) {
  print "$s->[0]\t$s->[1]\n";
}
