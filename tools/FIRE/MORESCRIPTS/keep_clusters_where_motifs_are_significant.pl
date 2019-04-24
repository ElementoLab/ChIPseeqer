BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $cnt = 0;
my %H   = ();
foreach my $r (@$a_ref) {
  my $m = shift @$r;
  foreach my $s (@$r) {
    if (!defined($H{$s})) {
      $H{ $s } = $cnt++;
    }
  }
}


$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref; print join("\t", @$r); print "\n";
foreach my $r (@$a_ref) {
  if (defined($H{$r->[1]})) {
    print "$r->[0]\t" . $H{$r->[1]} . "\n"; 
  }
}



