BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

  print Sets::myre2scanace($ARGV[0]);
  
  print "\n\n";



