open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if ($a[2] eq 'gene') {
    my $s = $a[8]; $s =~ /ID\=.+2\:(.+?)\;/;
    my $l = ($a[6] eq "+"?1:-1);
    print "$1\t$a[0]\t$a[3]\t$a[4]\t$l\n";
  }
  
}
close IN;
