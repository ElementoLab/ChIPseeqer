my %INDEX = ();

foreach my $f (@ARGV) {
  
  open IN, $f;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;

    if (!defined($INDEX{ $a[0] }{ $a[1] })) {
      print "$a[0]\t$a[1]\n";
      $INDEX{ $a[0] }{ $a[1] } = 1;
      $INDEX{ $a[1] }{ $a[0] } = 1;
    }

  }
  close IN;

}
