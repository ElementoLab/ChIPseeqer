open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $t = $a[1] + length($a[4]);
  print "$a[0]\t$a[1]\t$t\n";
}
close IN;

