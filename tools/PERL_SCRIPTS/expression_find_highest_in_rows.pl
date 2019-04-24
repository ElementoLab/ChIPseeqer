BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $cnt = 0;
foreach my $r (@$a_ref) {
  if ($cnt == 0) {
    print join("\t", @$r); print "\n";
  } else {
    
    my $g = shift @$r;
    
    my $n = scalar(@$r);
    
    my $max   = -10000000;
    my $max_i = undef;
    
    for (my $i=0; $i<$n; $i++) {
      if ($r->[$i] > $max) {
	$max_i = $i;
	$max   = $r->[$i];
      }
    }
    
    if (defined($max_i)) {
      print "$g";
      
      for (my $i=0; $i<$n; $i++) {
	if ($i == $max_i) {
	  print "\t1";
	}	else {
	  print "\t0";
	}
      } 
      
      print "\n";
    }
  }

  
  $cnt ++;
    
}

