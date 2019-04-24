#BEGIN{ $home = `echo \$HOME`; chomp $home}
#use lib "$home/PERL_MODULES";
use lib qw(/home/elemento/PERL_MODULES);

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_clu = $ta->getArray();


my $h_ref_zyg = Sets::getIndex($ARGV[2]);



$ta->loadFile($ARGV[1]);
my $a_ref_gen = $ta->getArray();
my $mm        = 5000;
foreach my $c (@$a_ref_clu) {

  print join("\t", @$c); 
  foreach my $g (@$a_ref_gen) {
    next if (!defined( $h_ref_zyg->{ $g->[0] } ));
    next if ($g->[1] ne $c->[0]);

    if (Sets::sequencesOverlap($c->[1], $c->[2], $g->[2], $g->[3])) {
      print "\tOV $g->[0]";
    } else {

      if ((abs($g->[2] - $c->[1]) < $mm) ||
	  (abs($g->[3] - $c->[1]) < $mm) || 
	  (abs($g->[2] - $c->[2]) < $mm) || 
	  (abs($g->[3] - $c->[2]) < $mm)) {
	
	if ($g->[4] == 1) {
	
	  if ($c->[2] < $g->[2]) {
	    print "\tUP $g->[0]";
	  } else {
	    print "\tDO $g->[0]";
	  }
	  
	} else {
	  if ($c->[1] > $g->[3]) {
	    print "\tDO $g->[0]";
	  } else {
	    print "\tUP $g->[0]";
	  }
	}
      }

    }
    
  }
  print "\n";

}

