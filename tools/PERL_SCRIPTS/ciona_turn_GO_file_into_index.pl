use lib qw(/home/elemento/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %G = ();
my %J = ();
foreach my $r (@$a_ref) {

  if ($r->[1] eq "Ciona intestinalis") {
    push @{ $J{ $r->[0] } }, $r->[2];
    $G{ $r->[2] } = $r->[3];
  }
}

foreach my $k (sort(keys(%J))) {
  print "$k\t" . join("\t", @{ $J{ $k } }); print "\n";
}

open OUT, ">GO-terms.txt";
foreach my $k (keys(%G)) {
  print OUT "$k\t$G{$k}\n";
}
close OUT;
