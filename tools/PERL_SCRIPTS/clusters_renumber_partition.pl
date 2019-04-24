BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;

print join("\t", @$r) . "\n";

$cnt ++;
foreach my $r (@$a_ref) {
  if (!defined($H{$r->[1]})) {
    $H{$r->[1]} = $cnt++;
  }
  
  $r->[1] = $H{ $r->[1] };
  print join("\t", @$r) . "\n";
}

