BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {
  $H{$r->[0]}{$r->[1]} ++;
}


foreach my $m (keys(%H)) {
  my $n = scalar(keys(%{$H{$m}}));
  my $f = $n / 16450;  
  print "$m\t$f\n";  
}
