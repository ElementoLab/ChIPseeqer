open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  next if ($l =~ /\#/);
  
  my @a = split /\t/, $l, -1;
  if ( ($a[1] eq "curated") && ($a[2] eq "Sequence")) {
    
    ($a[8]) = $a[8] =~ /Sequence \"(.+?)\"/;
    
    $a[6] = ($a[6] eq "+"?1:-1);

    print "$a[8]\t$a[0]\t$a[3]\t$a[4]\t$a[6]\n";
    
  }
  
}
close IN;
