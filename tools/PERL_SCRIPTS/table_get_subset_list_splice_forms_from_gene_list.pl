BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;

my $h_ref = Sets::getIndex($ARGV[0]);

my $ta    = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my @a = split /\-/, $r->[0];
  if (defined($h_ref->{$a[0]})) {
    print join("\t", @$r) . "\n";
  }
}


