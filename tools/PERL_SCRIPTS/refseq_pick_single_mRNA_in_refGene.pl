BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my $n = $r->[0];
  
  if ($r->[3] eq '+') {
    $st = $r->[4];
  } else {
    $st = $r->[5];
  }
  
  if (!defined( $H{ $r->[2] }{ $st }) ) {
    print join("\t", @$r); print "\n";
  }

  $H{ $r->[2] }{ $st } = 1;
}

