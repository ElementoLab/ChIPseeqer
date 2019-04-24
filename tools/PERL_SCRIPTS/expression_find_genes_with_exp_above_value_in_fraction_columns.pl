open IN, $ARGV[0] or die "Cannot open file";
my $l = <IN>;

print "GENE\tEXP\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $n = shift @a;
  
  my $num = 0;
  my $tot = 0;
  foreach my $r (@a) {
    if ($r >= $ARGV[1]) {
      $num++;
    }	
    $tot ++;
  }

  if ($num/$tot > $ARGV[2]) {
    print "$n\t1\n";
  } else {
    print "$n\t0\n";
  }


}
close IN;

