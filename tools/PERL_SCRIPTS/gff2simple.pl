open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  next if ($l =~ /^\#/);
  chomp $l;
  my @a = split /\t/, $l, -1;

  next if ($a[2] ne "CDS");
  
  my $chr = "H37Rv";
  my $ann = $a[8];

  my $tsname = undef;
  if ($ann =~ /^locus_tag=(.+?)\;/) {
    $tsname = $1;
  } elsif ($ann =~ /ID/) {
    my @b = split /\:/, $ann;
    my ($lo) = $ann =~ /locus_tag=(.+?)\;/;
    $tsname = "$lo-$b[1]";
  }
  
  if ($a[6] eq '+') {
    $a[6] = 1;
  } else {
    $a[6] = -1;
  }
  
  print "$tsname\t$chr\t$a[3]\t$a[4]\t$a[6]\n";

}
close IN;

