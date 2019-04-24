open IN, $ARGV[0];
my $row = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if ($row > 0) {
    $a[$ARGV[1]] = - $a[$ARGV[1]];
  }
  print join("\t", @a) . "\n";
  $row++;
}
close IN;
