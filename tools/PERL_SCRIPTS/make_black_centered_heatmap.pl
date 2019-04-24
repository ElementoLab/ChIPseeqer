my @c1 = split /\,/, $ARGV[0];
my @bl = (0,0,0);
my @c2 = split /\,/, $ARGV[1];


for (my $i=32; $i>=0; $i--) {
  
  my @b =();
  for (my $j=0; $j<3; $j++) {
    my $inc = ($c1[$j] - $bl[$j])/32;
    my $val = $i * $inc;
    push @b, sprintf("%3.2f", $val/255);
  }
  print join("\t", @b);
  print "\n";
  
  
}


for (my $i=1; $i<=32; $i++) {
  
  my @b =();
  for (my $j=0; $j<3; $j++) {
    my $inc = ($c2[$j] - $bl[$j])/32;
    my $val = $i * $inc;
    push @b, sprintf("%3.2f", $val/255);
  }
  print join("\t", @b);
  print "\n";
  
  
}
