open IN, $ARGV[0];
my $i = 0;

while (my $l = <IN>) {

  if ( ($i % 2) == 0) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    
    my @b = ();
    for(my $i=0; $i<@a; $i+=2) {
      push @b, $a[$i];
    }
    print join("\t", @b) . "\n";
  }
  $i++;
}
