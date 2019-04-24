my $g = shift @ARGV;


my %H = ();
foreach my $f (@ARGV) {
  open IN, $f;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    $H{$a[0]}{$f} = $a[1];
  }
close IN;


  
}


open IN, $g;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  print "$a[0]";
  foreach my $f (@ARGV) {
    print "\t" . $H{$a[0]}{$f};
  }
  print "\n";
}
close IN;
