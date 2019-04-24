open IN, $ARGV[0];
my $l = <IN>;
print $l;
my $cnt = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (!defined($H{$a[1]})) {
    $H{$a[1]} = $cnt ++;
  }

  print "$a[0]\t$H{$a[1]}\n";

}
close IN;

open OUT, ">$ARGV[0].idx.txt";
foreach my $k (sort(keys(%H))) {
  print OUT "$k\t$H{$k}\n";
} 
close OUT;
