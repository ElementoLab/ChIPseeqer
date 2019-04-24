BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

print Sets::count_file_lines($ARGV[0]);
print "\n";
