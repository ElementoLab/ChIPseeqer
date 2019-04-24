open IN, $ARGV[0];

my $i = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;

  for (my $j=0; $j<@a; $j++) {
    if ($a[$j] eq '') {
      print "Emptry entry at row $i, col $j\n";
    }
  }
  $i ++;
}
close IN;
