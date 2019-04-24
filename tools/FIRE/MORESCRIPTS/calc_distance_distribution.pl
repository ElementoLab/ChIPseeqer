use Sets;
use Table;

if (@ARGV == 0) {
  die "Usage : perl calc_distance_distribution.pl re profiles nbbins intl intr\n";
}

my $re   = $ARGV[0];
my $pr   = $ARGV[1];

my @POS  = (); 

open IN, $pr;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[0] eq $re) {
    push @POS, $a[2];
  }
}
close IN;

my $a = Sets::getSpecifiedDistribution(\@POS, $ARGV[2], $ARGV[3], $ARGV[4]);

foreach my $r (@$a) {
  print join("\t", @$r); print "\n";
}
