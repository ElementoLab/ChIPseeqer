open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  print ">$a[0]\n$a[1]\n\n";

}
close IN;

