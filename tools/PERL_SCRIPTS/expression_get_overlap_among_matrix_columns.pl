BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Hypergeom;

my $t = log(1.5)/log(2.0);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $n     = $ta->getNbColumns();

my @cols  = ();
my @labe  = ();
for (my $i=1; $i<$n; $i++) {
  my $a_ref = $ta->getColumn($i);
  my $l     = shift @$a_ref;
  push @labe, $l;
  push @cols, $a_ref;
}

for (my $i=0; $i<$n-1; $i++) {
  
  print $labe[$i];

  for (my $j=0; $j<$n-1; $j++) {

    my $s1 = 0;
    my $s2 = 0;
    my $ov = 0;
    my $nn = 0;
    for (my $k=0; $k<@{$cols[$i]}; $k++) {
      if ($cols[$i]->[$k] > $t) {
	$s1 ++;
      }
      if ($cols[$j]->[$k] > $t) {
	$s2 ++;
      }
      if (($cols[$i]->[$k] > $t) && ($cols[$j]->[$k] > $t)) { 
	$ov ++;
      }
      
      $nn++;
    }

    my $p = Hypergeom::cumhyper($ov, $s1, $s2, $nn);
    print sprintf("\t%1.1e", $p);
    
  }

  print "\n";

}

