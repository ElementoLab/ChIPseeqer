BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $w = 200;
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
 
    my $l = length($s);
    
    for (my $i=0; $i<$l-$w; $i++) {
      
      my $ss = substr($s, $i, $w);

      for (my $j=$i; $j<$i+$w; $j++) {
	
	
      
      }
      
    }
    
}
