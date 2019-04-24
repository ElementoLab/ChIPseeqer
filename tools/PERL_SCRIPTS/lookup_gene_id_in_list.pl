BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref1 = $ta->getArray();

$ta->loadFile($ARGV[1]);
my $a_ref2 = $ta->getArray();


foreach my $r1 (@$a_ref1) {

  my $t = 0;
  foreach my $r2 (@$a_ref2) {

    if ($r1->[0] eq $r2->[0]) {
      print "$r1->[0] found.\n"; 
      $t = 1;
      last;
    }

  }
  
  if ($t == 0) {
    print "$r1->[0] not-found\n";
  }
  
}

