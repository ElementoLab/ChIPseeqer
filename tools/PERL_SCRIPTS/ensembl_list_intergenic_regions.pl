BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

for (my $i=0; $i<@$a_ref; $i++) {
  my $r = $a_ref->[$i];

  # find the non-overlapping gene that's closest 5'
  # if there exists an overlapping gene, stops
  if ($r->[4] == 1) {
    for (my $j=0; $j<@$a_ref; $j++) {
      next if ($i == $j);
      
      if (Sets::sequenceOverlap()) {
	
	$dead = 1;
      }
      
    }
  }
  
  
  
}

