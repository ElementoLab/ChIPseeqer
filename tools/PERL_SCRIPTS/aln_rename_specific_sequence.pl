open IN, $ARGV[0];
while (my $l = <IN>) {
  #chomp $l;
  my $txt = $ARGV[1];

  if ($l =~ /(^$txt\ +)/) {
    my $t  = $1;
    my $tc = $t;
    my $l1 = length($t);
    $t =~ s/\ +$//;

    my $nt = $ARGV[2];
    my $l2 = length($nt);
    
    my $toadd = $l1 - $l2;
    $nt .= " " x $toadd;

    $l =~ s/$tc/$nt/;

    print $l;

  } else {
    print $l;
  }	

}
close IN;

