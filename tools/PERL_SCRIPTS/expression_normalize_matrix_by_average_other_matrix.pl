BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();


#
# calculate average
#
shift @$a_ref;
my @AVG = ();
my @NAM = ();
foreach my $r (@$a_ref) {
  my $g = shift @$r;
  push @NAM, $g;
  my $a = Sets::average($r);
  push @AVG, $a;
}


#
#  read in first matrix
#
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


#
# calculate average
#
my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

my $i = 0;
foreach my $r (@$a_ref) {
  
  my $n = shift @$r;
  die "toto\n" if ($NAM[$i] ne $n);
  print "$n";
  foreach my $s (@$r) {
    $s = sprintf("%3.2f", Sets::log2( $s / $AVG[$i] ));
    print "\t$s";
  }
  print "\n";

  $i++;
}

