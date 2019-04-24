while (my $l = <STDIN>) {
  chomp $l;
  my @a = split /\t/, $l, -1; 
  print "$a[0]\t$a[1]\t$a[2]-$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\n";

}
