open IN, $ARGV[0];
while (my $l = <IN>) 
{
  if ($l =~ /\# Name/) {
    my $l = <IN>;
    chomp $l;

    open OUT, ">$l.txt";
    #print "$l\n";
    
    my $l = <IN>;
    my $l = <IN>;
    for (my $i=0; $i<4; $i++) {
      my $l = <IN>; $l =~ s/\|//g;
      print OUT $l;
    }
    close OUT;
  }

}
