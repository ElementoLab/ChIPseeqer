
my %H = ();
my %G = ();
foreach my $f (@ARGV) {

  open IN, $f;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my $id = "$a[0]-$a[1]";
    $H{$id} = 1;
    push @{$G{$id}}, @a;
  }
  close IN;
  
  
}


foreach my $g (keys(%H)) {

  if ($H{$g} == @ARGV) {
    print join("\t", @{$G{$g}}) . "\n";
  }
  
}


