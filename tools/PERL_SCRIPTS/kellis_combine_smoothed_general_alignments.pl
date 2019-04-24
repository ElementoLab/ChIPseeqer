BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;

die "define w\n" if (!defined($ARGV[1]));
die "define conservation threshold\n" if (!defined($ARGV[2]));

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
    
    
    my @ah = split //, $lh;
    my @am = split //, $lm;
    my @ar = split //, $lr;
    my @ac = split //, $lc;


    my $n  = scalar(@ah);


    #my @ah_n = ();
    #my @am_n = ();
    #my @ar_n = ();
    #my @ac_n = ();
    # 
    #for (my $i=0; $i<$n; $i++) {
    #  if (($ah[$i] ne "-") && ($am[$i] ne "-")) {
    #	push @ah_n, $ah[$i];
    #	push @am_n, $am[$i];
    #  }
    #}
    
    
    #@ah = @ah_n;
    #@am = @am_n;

    $n  = scalar(@ah);


    # conservation vector
    my @as = ();
    for (my $i=0; $i<$n; $i++) {
      
      $as[$i] = 0;
      $as[$i] ++ if ($ah[$i] eq $am[$i]); 
      $as[$i] ++ if ($ah[$i] eq $ar[$i]); 
      $as[$i] ++ if ($ah[$i] eq $ac[$i]); 
    
    }


    my @asmo = ();    
    my $w = $ARGV[1];

    for (my $i=0; $i<$n; $i++) {

      if (($ah[$i] eq 'N') || ($ah[$i] eq '#') || ($ah[$i] eq '-')) {

	print "N"; 

      } else {
	
	my $cnt = 0;
	my $c_w = 0;

	for (my $j=$i-1; $j>=Sets::max($i-$w/2, 0); $j--) {
	  last if ($ah[$j] eq '#');
	  $cnt += $as[$j];
	  $c_w ++;
	
	}
	
	
	
	for (my $j=$i; $j<Sets::min($i+$w/2, $n); $j++) {
	  last if ($ah[$j] eq '#');
	  $cnt += $as[$j];
	  $c_w ++;
	
	}
      
	

	$asmo[ $i ] = $cnt / $c_w;
	
	if ($asmo[ $i ] >= $ARGV[2]) {
	  print "$ah[$i]";
	} else {
	  print "N";
	}
      }
    }
    
    print "\n";
        
  }
  
 

}
close IN;
