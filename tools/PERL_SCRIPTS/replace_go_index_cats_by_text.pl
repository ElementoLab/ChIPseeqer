BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;

my $h_ref = undef;
if ($ARGV[1]) {
  $ta->loadFile($ARGV[1]);
  $h_ref = $ta->getIndexKV(0,1);
}


$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my $n = shift @$r;
  if (defined($h_ref)) {
    foreach my $s (@$r) {
      $s = $h_ref->{$s};
    }
  }
  print "$n\t" . join("/", @$r) . "\n";
}

