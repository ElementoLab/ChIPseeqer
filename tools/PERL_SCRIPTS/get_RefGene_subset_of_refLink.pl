open IN, $ARGV[0];

my %H = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $H{$a[1]} = 1;
}
close IN;

open IN, $ARGV[1];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if (defined($H{$a[2]})) {
    print "$a[2]\t$a[3]\t$a[0]\t$a[1]\n";
  }
}
close IN;




