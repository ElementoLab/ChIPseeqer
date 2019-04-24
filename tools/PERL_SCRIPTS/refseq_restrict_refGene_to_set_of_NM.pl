BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;


my $h_ref = Sets::getIndex($ARGV[0]);



my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if (defined($h_ref->{ $r->[1] })) {
    print join("\t", @$r) . "\n";
  }
}

