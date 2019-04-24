while (my $l = <STDIN>) {
  
  chomp $l;

  my @a = split /\t/, $l, -1;

  $H{ $a[4] }{ $a[0] } ++;
  $C{ $a[4] } ++;
}


foreach my $k (keys(%C)) {
  print "$k\t$C{$k}\t";

  my @a  = (); 
  foreach my $kk (keys( %{ $H{$k} } )) {
    #print "$kk/$H{$k}{$kk} ";
    my @a_tmp = ($kk, $H{$k}{$kk});
    push @a, \@a_tmp;
  } 

  @a = sort { $b->[1] <=> $a->[1] } @a;

  foreach my $r (@a) {
    print "$r->[0]/$r->[1] ";
  }
  
  print "\n";

}
