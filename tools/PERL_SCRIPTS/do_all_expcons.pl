BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $d = "$ARGV[0]\_EXPCON";
mkdir $d if (! -e $d);

my $a_ref = $ta->getColumn(0);

my $cnt = 1;
foreach my $r (@$a_ref) {
  print "Processing $r\n";
  my $todo = "perl  ./expcons_read_fa_pp_all.pl chrom.txt \"$r\" $d/$cnt";
  system($todo);
  $cnt ++;
}
