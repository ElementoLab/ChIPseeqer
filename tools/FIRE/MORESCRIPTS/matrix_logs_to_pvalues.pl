BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {
  my $n = shift @$r;
  print "$n";
  foreach my $s (@$r) {
    my $l = undef;
    if ($s > 0) {
      $l = sprintf("%3.2e", exp( (- $s) * log(10.0)  ));
      $l = "$l(o)";
    } else {
      $l = sprintf("%3.2e", exp( (+ $s) * log(10.0)  ));
      $l = "$l(u)";
    }
    print "\t$l";
  }
  print "\n";
}

