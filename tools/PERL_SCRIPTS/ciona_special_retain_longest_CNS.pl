use lib qw(/home/elemento/PERL_MODULES);

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $n     = scalar(@$a_ref);

for (my $i=0; $i<$n; $i++) {
  
  next if ($TAKEN[$i] == 1);

  for (my $j=$i+1; $j<$n; $j++) {
    
    if (($a_ref->[$i]->[0] eq $a_ref->[$j]->[0]) && (Sets::sequencesOverlap($a_ref->[$i]->[3],
			       $a_ref->[$i]->[4],
			       $a_ref->[$j]->[3],
			       $a_ref->[$j]->[4]))) {
      $TAKEN[ $j ] = 1;
	
    }
  

  }
  

}

for (my $i=0; $i<$n; $i++) {
  if (!defined($TAKEN[$i])) {
    print join("\t", @{ $a_ref->[$i] } );
    print "\n";
  }
}

