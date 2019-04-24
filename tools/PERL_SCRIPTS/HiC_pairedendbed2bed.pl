open IN, $ARGV[0] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\ /, $l, -1;
  my @b = split /\,/, $a[10];

  my $r1_st = $a[1];
  my $r1_en = $a[1] + $b[0];

  my $r2_st = $a[2] - $b[1];
  my $r2_en = $a[2];

  print "$a[0]\t$r1_st\t$r1_en\tX\t0\t$a[5]\n";
  print "$a[0]\t$r2_st\t$r2_en\tX\t0\t$a[5]\n";

}
close IN;

