BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r = shift @$a_ref;

print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {
  my $g = shift @$r;
  print "$g";
  foreach my $s (@$r) {
    print "\t" . sprintf("%4.1f", $s);
  }
  print "\n";
}

