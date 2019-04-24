BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
#use lib qw(/home/elemento/PERL_MODULES);

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $cnt = 0;
foreach my $r (@$a_ref) {
  if ($cnt == 0) {
    print join("\t", @$r); print "\n";
  } else {
    my $g = shift @$r;
    my $a_ref_ranks = Sets::get_rank_from_array($r);
    print "$g\t" . join("\t", @$a_ref_ranks);
    print "\n";
  }
  $cnt ++;
}


