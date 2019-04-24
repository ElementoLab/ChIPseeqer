BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my %RANKS = ();
my $cnt   = 1;
foreach my $r (@$a_ref) {
  my $cc = Sets::getComplement($r->[0]);  
  $RANKS{ $r->[0] } = $cnt;
  $RANKS{ $cc     } = $cnt;
  $cnt++;
}


$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $cnt = 1;
foreach my $r (@$a_ref) {
  print "$r->[0]\t$cnt\t$RANKS{$r->[0]}\n"; 
  $cnt ++;
}


