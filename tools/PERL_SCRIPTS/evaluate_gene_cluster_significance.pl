BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";




my $l_total = $ARGV[0];
my $n_total = $ARGV[1];
my $l_win   = $ARGV[2];
my $n_win   = $ARGV[3];

srand;

for (my $i=0; $i<1000000; $i++) {
  
  my @POS = ();
  for (my $j=0; $j<$n_total; $j++) {
    my $p = 1 + int(rand($l_total));  #print "$p\n"; 
    push @POS, $p;
  }
  @POS = sort { $a <=> $b } @POS;

  #
  # find out if there is window of size $l_win
  #   containing at least $n_win occurrences 
  #  

  #print "$i\n";
  #print join("\t", @POS); print "\n";

  for (my $j=0; $j<$n_total-1; $j++) {
    my $p1 = $POS[ $j ];
    my $w1 = 1;
    my $out = 0;
    for (my $k=$j+1; $k<$n_total; $k++) {
      my $p2 = $POS[ $k ];
      if ( ($p2 - $p1) < $l_win ) { 
	if ( ($k-$j+1) >= $n_win ) {
	  print "WINDOW!\n";
	  $out = 1;
	  last;
	}  
      }
    }
    if ($out == 1) {
      last;
    }
  }

}
