BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref; print join("\t", @$r); print "\n";

my %H = ();
foreach my $r (@$a_ref) {
  $H{ $r->[0] } = 1;
}


$ta->loadFile($ARGV[1]);

my $r = shift @$a_ref; 

my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if (defined($H{ $r->[0] })) {
    print join("\t", @$r); print "\n";
  }
}

