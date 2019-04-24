BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

for (my $i=0; $i<@$a_ref; $i++) {

  open OUT, ">tmp";
  # output a file without this guy
  for (my $j=0; $j<@$a_ref; $j++) {
    print OUT join("\t", @$a_ref->[$i]) . "\n" if ($i != $j);
  }
  close OUT;
  
  my $todo ='perl ../PBS_mi_find.pl --files=tmp --fa=124_Dmel_Enc.fa --quantized=1 --target_dir=. --submit=0 --rna=0 --mbins_dist=4 --add=1';
  system("$todo");
}

