BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use strict;
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %O5 = ();
my %O3 = ();
foreach my $r (@$a_ref) {

  if ($r->[1] eq 'o5') {
    $O5{$r->[0]} = $r->[4];
  }
  
  if ($r->[1] eq 'o3') {
    $O3{$r->[0]} = $r->[4];
  }

}


foreach my $k (keys(%O5)) {
  print "$k\t" . ($O5{$k}-$O3{$k}) . "\n";
}
