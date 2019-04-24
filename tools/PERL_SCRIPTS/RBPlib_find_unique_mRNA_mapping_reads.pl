open IN, $ARGV[0] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  if (($a[1] =~ /mRNA/) && ($a[7] == 1)) {
    print "$l\n";
  }

}
close IN;

