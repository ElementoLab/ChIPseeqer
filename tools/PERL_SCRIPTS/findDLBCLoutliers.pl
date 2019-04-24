open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $r1 = ($a[3] + 1) / ($a[1] + $a[2] + 2);
  
  print join("\t", @a) . "\t$r1\n";

}
close IN;

