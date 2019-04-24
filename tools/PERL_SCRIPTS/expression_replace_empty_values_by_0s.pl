


my @a = split /\t/, $l, -1;

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  foreach my $r (@a) {
    if ($r eq "") {
      $r = 0;
    }
  }
  print join("\t", @a) . "\n";
}
close IN;

