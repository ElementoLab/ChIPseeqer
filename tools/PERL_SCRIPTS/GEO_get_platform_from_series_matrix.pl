open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  
  if ($a[0] =~ /\!Series_platform_id/) {
    #print "$l\n";
    $a[1] =~ s/\"//g;
    print "$a[1]\n";
    last;
  }

}
close IN;
