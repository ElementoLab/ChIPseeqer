open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  if ($l =~ /^\>/) {
    print "$l\n";
  } else {
    
    my $lh = <IN>;
    my $lm = <IN>;
    my $lr = <IN>;
    my $lc = <IN>;
    
    #my $ls = <IN>;
 
    my @ah = split //, $lh;
    my @am = split //, $lm;
    
    my $n  = scalar(@ah);

    for (my $i=0; $i<$n; $i++) {
      if ($ah[$i] eq $am[$i]) {
	if ($ah[$i] eq "#") {
	  print "N";
	} elsif ($ah[$i] ne "-") {
	  print "$ah[$i]";
	}
      } else {
	print "N";
      }	
    }
    
    print "\n";

    #<STDIN>;
  }

 

}
close IN;
