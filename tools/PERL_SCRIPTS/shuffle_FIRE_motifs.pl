BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {


  $r->[0] = Sets::shuffle_re($r->[0]);

  print join("\t", @$r) . "\n";

}


  
