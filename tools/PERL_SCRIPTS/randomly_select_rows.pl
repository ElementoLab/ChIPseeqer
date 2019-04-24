BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r); print "\n";
foreach my $r (@$a_ref) {
  
  my $v = rand();
  if ($v <= 0.2) {
    print join("\t", @$r); print "\n";
  }
  
}

