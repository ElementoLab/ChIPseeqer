open IN, $ARGV[1];
my %H = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $H{$a[0]} = 1;
}
close IN;


open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (!defined($H{$a[0]})) {
    print "$a[0]\n";
  }

}
close IN;


