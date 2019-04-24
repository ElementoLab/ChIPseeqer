BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {
  my $n = shift @$r;
  foreach my $s (@$r) {
    push @{ $H{ $s } }, $n if (!Sets::in_array($n, @{ $H{ $s } }));
  }
}


foreach my $r (keys(%H)) {
  print "$r\t" . join("\t", @{ $H{ $r } }) . "\n";
}
