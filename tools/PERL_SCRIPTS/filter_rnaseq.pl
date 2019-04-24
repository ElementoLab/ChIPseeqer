open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (($a[2] < $a[1]) && ($a[3] < $a[1])) {
    print "$l\n";
    #print "$a[0]\t1\n";
  } else {
    #print "$a[0]\t0\n";
  }

}
close IN;

