open IN, $ARGV[0];
my $i = 0;
my $j = 0;
while (my $l = <IN>) {

  if ( ($i % 2) == 0) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    print "$j\t$a[1]\n";
    $j++;
  }
  $i++;
}
