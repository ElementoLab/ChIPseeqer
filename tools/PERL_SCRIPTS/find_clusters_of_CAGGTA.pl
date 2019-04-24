BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $n     = scalar( @$a_ref );

my $i = 0;
while ($i<$n) {  

  my $p1 = $a_ref->[ $i ]->[1];
  #print $a_ref->[$i]->[0] . "\t" . $a_ref->[$i]->[1] . "\n";
  my $m  = 1;
  my $c1 = $a_ref->[ $i ]-> [0];
  my $l_cur = 1;

  my $j  = $i+1;
  while ( $j < $n ) {
    my $p2 = $a_ref->[ $j ]->[1];    
    my $c2 = $a_ref->[ $j ]-> [0];
    
    last if ($c1 ne $c2);

    my $l  = $p2 - $p1;
    if ($l < 1500) {
      $m ++;
      $l_cur = $l;
      $p2_cur = $p2;
      #print "$l\t$m\n";
    } else {      
      if (($m > 1) && ($m/$l_cur >= 2/1000) && !Sets::sequencesOverlap($p1, $p2, $p1_last, $p2_last)) {

	
	print "$c1\t$p1\t$p2_cur\t$l_cur\t$m\n";
	$p1_last = $p1;
	$p2_last = $p2;
      }
      last;
    }
    $j ++;
  } 
  $i ++;
}
