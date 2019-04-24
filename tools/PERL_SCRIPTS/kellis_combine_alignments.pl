open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  
  next if ($l eq "");

  

  if ($l =~ /^\>/) {

    print "$l\n";
    
    my $lh = <IN>;
    my $lm = <IN>;
    my $lr = <IN>;
    my $lc = <IN>;
    
    #my $ls = <IN>;
 
    my @ah = split //, $lh;
    my @am = split //, $lm;
    my @ar = split //, $lr;
    my @ac = split //, $lc;
    
    my $n  = scalar(@ah);

    for (my $i=0; $i<$n; $i++) {
      if (($ah[$i] eq $am[$i]) && ($ah[$i] eq $ar[$i]) && ($ah[$i] eq $ac[$i])) {
	if ($ah[$i] eq "#") {
	  print "N";
	} elsif ($ah[$i] eq "-") {
	  print "N";
	} else {
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
