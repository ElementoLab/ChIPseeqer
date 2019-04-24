open IN,$ARGV[0];
while (my $l = <IN>) {
  chomp $l;

  if ($l =~ /setrgbcolor$/) {
    my @a = split /\ /, $l, -1;
    if ($a[0] != int($a[0])) {
      $a[0] = int(0.5+$a[0]*100) / 100;      
    }
    if ($a[1] != int($a[1])) {
      $a[1] = int(0.5+$a[1]*100) / 100;      
    }
    if ($a[2] != int($a[2])) {
      $a[2] = int(0.5+$a[2]*100) / 100;      
    }
    print join(" ", @a) . "\n";
  } else {
    print "$l\n";
  }


}
close IN;
