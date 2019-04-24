BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  
  next if ($l eq "");

  

  if ($l =~ /^\>/) {

    print "$l\n";
    
    my $lh = <IN>; $lh = uc($lh);
    my $lm = <IN>; $lm = uc($lm);
    my $lr = <IN>; $lr = uc($lr);
    my $lc = <IN>; $lc = uc($lc);
    
    #my $ls = <IN>;
 
    my @ah = split //, $lh;
    my @am = split //, $lm;
    my @ar = split //, $lr;
    my @ac = split //, $lc;
    
    my $n  = scalar(@ah);
    
    my @as = ();

    for (my $i=0; $i<$n; $i++) {
      
      $as[$i] = 0;
      $as[$i] ++ if ($ah[$i] eq $am[$i]); 
      $as[$i] ++ if ($ah[$i] eq $ar[$i]); 
      $as[$i] ++ if ($ah[$i] eq $ac[$i]); 
    
    }

    my @asmo = ();
    
    my $w = 10;
    for (my $i=0; $i<$n; $i++) {
      my $cnt = 0;
      for (my $j=Sets::max($i-$w/2, 0); $j<Sets::min($i+$w/2, $n); $j++) {
	$cnt += $as[$j];
      }


      $asmo[ $i ] = $cnt / $w;

      #my $c = int( 0.5 + $asmo[$i]);
      #print "$c$ah[$i]";
		   
      if ($asmo[ $i ] >= 2.5) {
	
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
