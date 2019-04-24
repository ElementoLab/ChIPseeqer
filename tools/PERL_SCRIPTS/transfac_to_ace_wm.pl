BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

open IN, $ARGV[0];

my @W_RAW = ();
my @nt = ('A', 'C', 'G', 'T');
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if ($l =~ /^P0/) {
    
    my $cntp = 0;
    while (my $m = <IN>) {
      chomp $m;
      last if ($m =~ /^XX/);
      
      my @b = split /\ +/, $m, -1;
      shift @b;
      pop @b;
      
      
      my $s = Sets::arraySum(\@b);

      my $mul = 100 / $s ;
      
      #print "mul=$mul\n";
      my $p = "";
      for(my $i=0; $i<@b; $i++) {
	$p .= $nt[$i] x int(0.5 + $b[$i]*$mul);
      }

      my @d = split //, $p;
      push @W_RAW, \@d;
      
      $cntp ++;
    }
    
    print "Motif 1\n";
    my $l = scalar(@{ $W_RAW[0] });
    for (my $i=0; $i<$l; $i++) {
      for (my $j=0; $j<$cntp; $j++) {
	print (defined($W_RAW[$j][$i])?$W_RAW[$j][$i]:'N');
      }
      print "\n";
    }
    print "*" x $cntp; print "\n";
    last;
  }

}
close IN;


  

