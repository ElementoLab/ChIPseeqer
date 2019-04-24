BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my %CAT = ();
foreach my $r (@$a_ref) {
  my $g = shift @$r;
  foreach my $c (@$r) {
    push @{ $CAT{ $c } }, $g;
  }
}

$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0,1);

foreach my $c (keys(%CAT)) {
  if (@{$CAT{$c}} <= 300) {
    print "$c\t" . $h_ref->{$c} . "\t" . join("\t", @{$CAT{$c}}) . "\n";
  }
}
