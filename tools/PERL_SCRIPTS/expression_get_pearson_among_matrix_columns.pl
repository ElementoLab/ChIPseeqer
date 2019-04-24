BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

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
    
    my $p = Sets::spearman($cols[$i], $cols[$j]);
    print sprintf("\t%4.3f", $p);
    
  }

  print "\n";

}

