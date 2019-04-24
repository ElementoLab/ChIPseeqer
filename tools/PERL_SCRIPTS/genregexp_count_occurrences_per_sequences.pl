BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %POS = ();
foreach my $r (@$a_ref) {
  push @{ $POS{ $r->[0] } }, $r->[2] if (!Sets::in_array($r->[2], @{ $POS{ $r->[0] } }));
}


foreach my $g (keys(%POS)) {
  my $n = scalar(@{ $POS{ $g } });
  print "$g\t$n\t" . join("\t", @{ $POS{ $g } }); print "\n";
}

