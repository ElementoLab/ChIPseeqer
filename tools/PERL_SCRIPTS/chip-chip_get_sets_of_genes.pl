BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $d = $ARGV[1];

my $n = $ta->getNbColumns();

my $a_ref_genes = $ta->getColumn(0);
shift @$a_ref_genes;

for (my $i=1; $i<$n; $i++) {
  my $a_ref_col = $ta->getColumn($i);
  
  my $tf = shift @$a_ref_col;
  my @GENES = ();
  for (my $j=0; $j<@$a_ref_col; $j++) {
    next if ($a_ref_col->[$j] eq "NaN");
    if ($a_ref_col->[$j] < 1e-2) {
      push @GENES, $a_ref_genes->[$j];
    }
  }
  
  $tf =~ s/\ /\_/g;
  
  Sets::writeSet(\@GENES, "$d/$tf.txt");

}


