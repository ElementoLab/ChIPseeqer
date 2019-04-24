BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;
use DB_File;

my $FILE = "$ARGV[0]";
my $INDEX = "$FILE.idx";
open INF, $FILE;
open IND, $INDEX;

my $row = Sets::line_with_index(*INF, *IND, $ARGV[1] - 1);

print "$row\n";
