BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $H = ();

my $r = shift @$a_ref;
print join("\t", @$r); print "\n";
foreach my $r (@$a_ref) {
  if ($H{$r->[1]} < 70) {
    if ($H{$r->[1]} >= 60) {
      print join("\t", @$r); print "\n";
    }
    $H{$r->[1]} ++;
  }

}

