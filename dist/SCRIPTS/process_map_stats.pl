open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /[\t\ ]+/, $l;


  my $r = $a[2]/$a[4];

    print "$a[0]\t$a[2]/$a[4]\t$r\n";


}
close IN;

